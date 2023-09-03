// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
SEXP MRR(Eigen::MatrixXd Y,
          Eigen::MatrixXd X,
          int maxit = 500,
          double tol = 10e-9,
          int cores = 1,
          bool TH = false,
          double NonLinearFactor = 0.0,
          bool InnerGS = false,
          bool NoInversion = false,
          bool HCS = false,
          bool XFA = false,
          int NumXFA = 2,
          double prior_R2 = 0.5,
          double gc_prior_df = 0.5, 
          double var_prior_df = 0.0, 
          double weight_prior_h2 = 0.0,
          double weight_prior_gc = 0.0,
          double PenCor = 0.0,
          double MinCor = 1.0,
          double uncorH2below = 0.0,
          double roundGCupFrom = 1.0,
          double roundGCupTo = 1.0,
          double roundGCdownFrom = 1.0,
          double roundGCdownTo = 0.0,
          double bucketGCfrom = 1.0,
          double bucketGCto = 1.0,
          double DeflateMax = 0.9,
          double DeflateBy = 0.005,
          bool OneVarB = false,
          bool OneVarE = false,
          bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXd xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXd ve = vy * (1-prior_R2);
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  Eigen::VectorXd vbInit = ((vy*prior_R2).array()/MSx.array());
  Eigen::VectorXd veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  double tmp;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      if(i>j){
        tmp = gc_prior_df * sqrt(vb(i,i)*vb(j,j));
        vb(i,j) = tmp;
        vb(j,i) = tmp;
      }
    }
  }
  
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXSX(k);
  Eigen::MatrixXd Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXd Sb = vb*var_prior_df;
  Eigen::VectorXd Se = ve*var_prior_df;
  Eigen::VectorXd iNp = (n.array()+var_prior_df-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  double bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0;
  int numit = 0;
  double logtol = log10(tol);
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Non-Linear weights for marker effects
  bool NonLinear = NonLinearFactor!=0.0;
  Eigen::MatrixXd W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXd iVeWj = iVe*1.0;
  Eigen::VectorXd tmpW(p);
  double maxW, minW;
  
  // Objects for other variance structures
  double gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(vb);
  Eigen::MatrixXd UDU(k,k);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      
      // System of equations - Traditional vs Stranden and Garrick 2009
      if(NoInversion){
        LHS = vb * (XX.row(J).transpose().array() * iVeWj.array()).matrix().asDiagonal(); 
        for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = (vb * (RHS.array() * iVeWj.array()).matrix()).array();
      }else{
        LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()).matrix();
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = RHS.array() * iVeWj.array();
      }
      
      // Update coefficient
      b0 = b.row(J)*1.0;
      for(int i=0; i<k; i++){ iVeWj(i) = iVe(i)*W(J,i); }
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()   ).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVeWj.array();
      
      // Inner GS
      if(InnerGS){
        b1 = b.row(J)*1.0;
        for(int i=0; i<k; i++){
          ri = InnerRGSvec[i];
          b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      }else{
        b1 = LHS.llt().solve(RHS); 
      }
      
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Update marker weights
    if(NonLinear){
      W = b.cwiseAbs();
      for(int j=0; j<k; j++){
        maxW = W.col(j).maxCoeff(); minW = W.col(j).minCoeff();
        tmpW = NonLinearFactor * (W.col(j).array()-minW)/(maxW-minW) + (1.0-NonLinearFactor);
        tmpW = tmpW.array() + (1.0-tmpW.mean());
        W.col(j) = tmpW.array();
      }
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    h2 = 1 - ve.array()/vy.array();
    // Proportion-based prior
    if(weight_prior_h2>0){for(int i=0; i<k; i++){gs = ve(i)*(1-weight_prior_h2) + weight_prior_h2*veInit(i); ve(i) = gs*1.0;}}
    // Single variance
    if(OneVarE){tmp = ve.array().mean(); for(int i=0; i<k; i++) ve(i) = tmp*1.0;}
    iVe = ve.array().inverse();
    iVeWj = iVe*1.0;
    
    //Genetic variance
    
    // Get tilde-hat
    if(TH){
      for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));
      }
      TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    }else{
      TildeHat = b.transpose()*tilde;
    }
    
    // Estimate variances and covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          if(TH){
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXSX(i)+var_prior_df);
          }else{
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrXSX(i)+var_prior_df);
          }
        }else{ // Covariances
          if(TH){
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrDinvXSX(i)+TrDinvXSX(j)+var_prior_df);
          }else{
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrXSX(i)+TrXSX(j)+var_prior_df);
          }
        }}}
    
    // TO-DO COMPUTE CORR HERE, AND ONLY ONCE !!!
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)= (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc_prior_df*weight_prior_gc;}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){ vb(i,j) =  GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
    
    // Heterogeneous Compound Symmetry
    if(HCS){
      gs = 0.0;
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i>j){gs += GC(i,j);}}}
      gs = gs/((k*(k-1))/2);
      // Extended Factor Analytics
    }else if(XFA){
      es.compute(GC);
      UDU = es.eigenvalues()[k] * es.eigenvectors().col(k) * es.eigenvectors().col(k).transpose();
      for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i] * es.eigenvectors().col(k-i) * es.eigenvectors().col(k-i).transpose();
      GC = UDU * 1.0; for(int i=0; i<k; i++){ GC(i,i)=1.0; };
    }
    
    // Monkeying with the correlations
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i!=j){
            // Zero'ing  Correlations
            if(MinCor<1.0){ if(GC(i,j)<MinCor){ GC(i,j) = 0.0; }}
            // Penalize Correlations
            if(PenCor>0.0){  GC(i,j) = tanh(PenCor*abs(GC(i,j)))*GC(i,j);} 
            // Round Down
            if(roundGCdownFrom<1.0){ if(GC(i,j)<roundGCdownFrom){ GC(i,j) = roundGCdownTo*1.0; }}
            // Round Up
            if(roundGCupFrom<1.0){ if(GC(i,j)>roundGCupFrom){ GC(i,j) = roundGCupTo*1.0; }}
            // Bucket round
            if(bucketGCfrom<1.0){ if(GC(i,j)>bucketGCfrom && GC(i,j)<bucketGCto  ){ GC(i,j) =  bucketMean*1.0; }}
            // Min H2
            if(uncorH2below>0.0){ if(h2(i)<uncorH2below || h2(j)<uncorH2below  ){ GC(i,j) = 0.0; }}
      }}}
    
    // RECONSTRUCT COVARIANCE HERE AND ONLY ONCE
    // Single variance of beta
    if(OneVarB){tmp = TildeHat.diagonal().mean(); vb = GC * tmp;  }else{
    // Regular covariance reconstruction
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){ vb(i,j) =  GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}
    
    // Bending
    if( !NoInversion || TH ){
      A = vb*1.0;
      for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= Deflate;}}}
      if(A.llt().info()==Eigen::NumericalIssue && Deflate>DeflateMax){
        Deflate -= DeflateBy; if(verbose) Rcpp::Rcout <<  Deflate;
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*Deflate;}}}}
      EVDofA.compute(A);
      MinDVb = EVDofA.eigenvalues().minCoeff();
      if( MinDVb < 0.0 ){ 
        if(verbose) Rcpp::Rcout << ".";
        inflate = inflate - MinDVb*1.001;
        A.diagonal().array() += inflate;}
      iG = A.inverse();
    }
    
    // Compute convergence and print status
    
    //cnv = log10((beta0.array()-b.array()).square().sum());
    cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
    CNV1(numit) = cnv; if(std::isnan(cnv)){ if(verbose){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n";} break;}
    CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
    CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
    
    // Print
    ++numit;
    if( verbose && numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){ if(verbose){Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n";} break;}
    if( numit == maxit && verbose){ Rcpp::Rcout << "Model did not converge\n"; }
    
  }
  
  // Fitting the model
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}

  // Resize convergence vectors
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){
    CNV1b(i)=CNV1(i);
    CNV2b(i)=CNV2(i);
    CNV3b(i)=CNV3(i);
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("vb")=vb,
                            Rcpp::Named("ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("bend")=Deflate,
                            Rcpp::Named("cnvB")=CNV1b,
                            Rcpp::Named("cnvH2")=CNV2b,
                            Rcpp::Named("cnvV")=CNV3b,
                            Rcpp::Named("b_Weights")=W,
                            Rcpp::Named("Its")=numit);
}

