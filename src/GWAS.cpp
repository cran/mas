// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
SEXP GWAS(Eigen::MatrixXd Y, Eigen::MatrixXd GEN, Eigen::MatrixXd M,
              int maxit = 500,double logtol=-8, int cores = 1, bool verb = true){
  
  // 
  // Null model - Multivariate SNP-BLUP
  // 
  
  // Start setup
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = GEN.cols(), f = M.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  if(verb){ Rcpp::Rcout << "Running absobtion of phenotypes (SY)\n"; }
  Eigen::MatrixXd y(n0,k), WM(n0,f), MU(f,k), iMM(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WM.col(j)=M.col(j).array()*Z.col(i).array();}
    iMM = (WM.transpose()*WM).inverse();
    MU.col(i) = iMM * WM.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WM*MU.col(i) ).array()*Z.col(i).array();
  }
  
  // Compute SX
  if(verb){  Rcpp::Rcout << "Running absobtion of marker scores (SZ)\n"; }
  iMM = (M.transpose()*M).inverse();
  for(int j=0; j<m; j++){ GEN.col(j) = (GEN.col(j) - M*(iMM*M.transpose()*GEN.col(j))).array(); }
  
  // Single value decomposition
  if(verb){ Rcpp::Rcout << "SVD of marker score matrix\n"; }
  Eigen::BDCSVD<Eigen::MatrixXd> svd(GEN, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::VectorXd D = svd.singularValues();
  Eigen::MatrixXd X = svd.matrixU() * D.array().matrix().asDiagonal();
  int p = X.cols();
  
  // Sum of squares of X
  if(verb){ Rcpp::Rcout << "Set starting values for coefficients and variances\n"; }
  Eigen::MatrixXd XX(p,k); 
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  Eigen::VectorXd MSX = XX.colwise().sum().array(); MSX=MSX.array()*iN.array();
  // Variances
  iN = (n.array()-f).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSX.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXX(k);
  Eigen::MatrixXd Dinv(p,k);
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  // Bending
  double deflate = 1.0, deflateMax = 0.9;
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0;
  // Convergence control
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0;
  int numit = 0;
  // Prior for stability
  double df0 = 1.1;
  Eigen::MatrixXd Sb = vb*df0;
  Eigen::VectorXd Se = ve*df0;
  Eigen::VectorXd iNp = (n.array()+df0-f).inverse();
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Loop
  if(verb){ Rcpp::Rcout << "Starting Gauss-Seidel loop\n"; }
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    //for(int J=0; J<p; J++){
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVe.array();
      if(k==1){ b1 = RHS.array()/LHS.array(); }else{
        b1 = LHS.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(RHS);}
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    //ve = ve.array() / (n.array()-f);
    ve = (ve.array()+Se.array()) * iNp.array();
    iVe = ve.array().inverse();
    h2 = 1 - ve.array()/vy.array();
    
    // Get tilde-hat
    for(int i=0; i<k; i++){
      Dinv.col(i) = (XX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
      TrDinvXX(i)  = (XX.col(i).transpose() * Dinv.col(i));}
    TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          //vb(i,i) = TildeHat(i,i)/TrDinvXX(i);
          vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXX(i)+df0);
        }else{ // Covariances
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrDinvXX(i)+TrDinvXX(j));
        }}}
    
    // Bending
    A = vb*1.0;
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= deflate;}}}
    if(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){  deflate -= 0.005;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}}
    EVDofA.compute(A);
    MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ 
      inflate = inflate - MinDVb*1.001;
      A.diagonal().array() += inflate;}
    iG = A.inverse();
    
    // Covariances
    cnv = log10((beta0.array()-b.array()).square().sum());  
    CNV1(numit) = cnv;
    CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
    CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
    if( std::isnan(cnv) && verb ){Rcpp::Rcout << "Numerical issue!!\n"; break;}
    
    // Print
    ++numit;
    if( numit % 100 == 0 && verb){ Rcpp::Rcout << "Iter: "<< numit << " || log10 Conv: "<< cnv << "\n"; } 
    if( numit > 10 && cnv<logtol && verb ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
    if( numit == maxit && verb ){ Rcpp::Rcout << "Model did not converge\n"; }
  }
  
  if(verb){ Rcpp::Rcout << "Fitting final model\n"; }
  // Fitting the model
  Eigen::MatrixXd hat = X*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( M * MU.col(i) + hat.col(i)).array(); }
  Eigen::MatrixXd beta = V*b;
  
  // Correlations
  if(verb){ Rcpp::Rcout << "Estimating correlations\n"; }
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  if(verb){  Rcpp::Rcout << "Convergence statistics\n"; }
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){ CNV1b(i)=CNV1(i);CNV2b(i)=CNV2(i);CNV3b(i)=CNV3(i);}
  
  // Null model Output
  Rcpp::List NullModelOutput = Rcpp::List::create(Rcpp::Named("FxdEffCoef")=MU,
                                                  Rcpp::Named("MarkerEffects")=beta,
                                                  Rcpp::Named("FittedValues")=hat,
                                                  Rcpp::Named("Heritability")=h2,
                                                  Rcpp::Named("GenCorrelations")=GC,
                                                  Rcpp::Named("VarBeta")=vb,
                                                  Rcpp::Named("VarResiduals")=ve,
                                                  Rcpp::Named("ConvergenceBeta")=CNV1b,
                                                  Rcpp::Named("ConvergenceH2")=CNV2b,
                                                  Rcpp::Named("ConvergenceVar")=CNV3b,
                                                  Rcpp::Named("NumOfIterations")=numit);
  
  // 
  // Membership GWAS
  // 
  if(verb){ Rcpp::Rcout << "Starting Membership GWAS\n"; }
  
  // Membership matrix and other things
  int df;
  double tmp, VarE, dff, MSB;
  Eigen::MatrixXd ZM(n0,f), PZM(n0,f), ZPZ(f,f);
  Eigen::MatrixXd QtlVar(m,k), QtlHer(m,k), Gamma(m,f), Wald(m,k), PopWald(m,f);
  Eigen::VectorXd mum(f), ZPy(f), gamma(f), iC(n0), tmpCoef(n0);
  Eigen::VectorXd yZg(n0), PyZg(n0);
  Eigen::MatrixXd MSEx(m,k);
  
  // Loop across traits
  for(int trait=0; trait<k; trait++){
    
    // Don't run if the base model broke
    if( std::isnan(cnv) && verb ){Rcpp::Rcout << "GWAS CANNOT RUN!!\n"; break;}
    
    
    // degrees of freedom
    df = n(trait)-1.0; 
    dff = 1.0/(df-f);
    
    // Bypassing P matrix
    iC = 1.0 / ( D.array().square() + (ve(trait)/vb(trait,trait)) ).array();
    
    // Genome screening
    if(verb){ Rcpp::Rcout << "Genome-wide screening (Trait "<< (trait+1) << ")\n"; }
    for(int J=0; J<m; J++){
      
      // Print marker  
      if( J>0 && J % 1000 == 0 && verb ){ Rcpp::Rcout << "Screening marker: "<< J << "\n"; } 
      
      // Create Shizhong's Z matrix
      for(int i=0; i<f; i++){
        ZM.col(i) =  M.col(i).array() * GEN.col(J).array();
        tmpCoef = iC.array().matrix().asDiagonal() * X.transpose() * ZM.col(i);
        PZM.col(i) = ZM.col(i).array() - (X*tmpCoef).array();}
      
      // Compute coefficient
      ZPZ = ZM.transpose() * PZM; 
      ZPZ.diagonal().array() += 0.1; // Added for algorithmic stability
      ZPy = ZM.transpose() * e.col(trait);
      gamma = ZPZ.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(ZPy);
      Gamma.row(J) = gamma.array();
      
      // Variance components and heritability
      yZg = (y.col(trait).array()-(ZM*gamma).array())* Z.col(trait).array();
      PyZg = e.col(trait).array()-(PZM*gamma).array();
      VarE = (yZg.transpose()*PyZg);
      VarE = VarE*dff;
      tmp = (gamma.transpose() * ZPZ * gamma).trace();
      MSB = tmp/f;
      
      // Store main var components
      MSEx(J,trait) = VarE;
      if(tmp>0){
        Wald(J,trait) = tmp/VarE;
        if((MSB-VarE)>0){
          QtlVar(J,trait) = f*(MSB-VarE)/df;
          QtlHer(J,trait) = QtlVar(J,trait)/(QtlVar(J,trait)+VarE);
        }else{
          QtlVar(J,trait)=0.0;
          QtlHer(J,trait)=0.0;
        }
      }else{
        QtlVar(J,trait)=0.0;
        QtlHer(J,trait)=0.0;
        Wald(J,trait)=0.0;
      }
      PopWald.row(J) = gamma.array().square() * ZPZ.diagonal().array()/VarE;}   
    
    // End GWAS loop
  }
  
  // GWAS output
  Rcpp::List GwasModelOutput = Rcpp::List::create(Rcpp::Named("Wald")=Wald,
                                                  Rcpp::Named("WaldByPop")=PopWald,
                                                  Rcpp::Named("QTLvar")=QtlVar,
                                                  Rcpp::Named("QTLh2")=QtlHer,
                                                  Rcpp::Named("Effect")=Gamma,
                                                  Rcpp::Named("Membership")=M,
                                                  Rcpp::Named("MSE")=MSEx,
                                                  Rcpp::Named("N")=n);
  // Final output
  Rcpp::List OutputList = Rcpp::List::create(Rcpp::Named("GBLUP")=NullModelOutput, 
                                             Rcpp::Named("GWAS")=GwasModelOutput);
  OutputList.attr("class") = "GenModel";
  return OutputList;
  
}


// [[Rcpp::export]]
SEXP MLM(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd Z,
         int maxit = 500, double logtol = -8, int cores = 1){
  
  // Basic info
  double df0 = 1.1; 
  if(cores!=1) Eigen::setNbThreads(cores);
  int k = Y.cols(), n0 = Y.rows(), f = X.cols(), p = Z.cols();
  
  // Incidence matrix W
  Eigen::MatrixXd W(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        W(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ W(i,j) = 1.0;}}}
  Eigen::VectorXd n = W.colwise().sum();
  Eigen::VectorXd iN = (n.array()-f).inverse();
  
  // Compute SY, SZ
  Eigen::MatrixXd y(n0,k), WX(n0,f), MU(f,k), iXX(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WX.col(j)=X.col(j).array()*W.col(i).array();}
    iXX = (WX.transpose()*WX).inverse();
    MU.col(i) = iXX * WX.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WX*MU.col(i) ).array()*W.col(i).array(); }
  iXX = (X.transpose()*X).inverse();
  for(int j=0; j<p; j++){ Z.col(j) = (Z.col(j) - X*(iXX*X.transpose()*Z.col(j))).array(); }
  
  // Sum of squares of Z
  Eigen::MatrixXd ZZ(p,k); 
  for(int i=0; i<p; i++){ ZZ.row(i) = Z.col(i).array().square().matrix().transpose() * W;}
  Eigen::VectorXd TrZSZ = ZZ.colwise().sum().array();
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // Variances
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/ (TrZSZ.array()*iN.array())  ).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd tilde = Z.transpose() * y;
  
  // Bending
  double deflate = 1.0, deflateMax = 0.9;
  Eigen::MatrixXd A = vb*1.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0;
  
  // Prior for stability
  Eigen::MatrixXd Sb = vb*df0;
  Eigen::VectorXd Se = ve*df0;
  Eigen::VectorXd iNp = (n.array()+df0-f).inverse();
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  Eigen::MatrixXd beta0(p,k);
  double cnv = 10.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (ZZ.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (Z.col(J).transpose()*e).array() + ZZ.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      
      // Update residuals
      e = (e-(Z.col(J)*(b1-b0).transpose()).cwiseProduct(W)).matrix();}
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    //ve = ve.array() * iN.array();
    ve = (ve.array()+Se.array()) * iNp.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){
        vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrZSZ(i)+df0);            }else{
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrZSZ(i)+TrZSZ(j)); }}}
    
    // Bending
    A = vb*1.0;
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= deflate;}}}
    if(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){
      deflate -= 0.005; Rcpp::Rcout <<  deflate;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}}
    EVDofA.compute(A);
    MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ 
      Rcpp::Rcout << ".";
      inflate = inflate - MinDVb*1.001;
      A.diagonal().array() += inflate;}
    iG = A.inverse();
    
    // Print status
    ++numit;
    cnv = log10((beta0.array()-b.array()).square().sum());
    if( std::isnan(cnv) ){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n"; break;}
    if( numit % 100 == 0 ){ Rcpp::Rcout << "Iter: "<< numit << " || log10 Conv: "<< cnv << "\n"; } 
    if(  cnv<logtol ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
    if( numit == maxit ){ Rcpp::Rcout << "Model did not converge\n"; }
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd hat = Z*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( X * MU.col(i) + hat.col(i)).array(); }
  
  // Genetic correlations
  Eigen::MatrixXd GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Name and create outputs
  Rcpp::List OutputList = Rcpp::List::create(Rcpp::Named("b")=MU,Rcpp::Named("u")=b,
                                             Rcpp::Named("hat")=hat,Rcpp::Named("h2")=h2,
                                             Rcpp::Named("GC")=GC,Rcpp::Named("bend")=deflate,
                                             Rcpp::Named("vb")=vb,Rcpp::Named("ve")=ve,
                                             Rcpp::Named("cnv")=cnv,Rcpp::Named("its")=numit);
  
  // Output
  OutputList.attr("class") = "PEGSmodel";
  return OutputList;
  
}
