// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
SEXP GWAS(Eigen::MatrixXf Y, Eigen::MatrixXf GEN, Eigen::MatrixXf M,
              int maxit = 500,float logtol=-8, int cores = 1, bool verb = true){
  
  // 
  // Null model - Multivariate SNP-BLUP
  // 
  
  // Start setup
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = GEN.cols(), f = M.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  if(verb){ Rcpp::Rcout << "Running absobtion of phenotypes (SY)\n"; }
  Eigen::MatrixXf y(n0,k), WM(n0,f), MU(f,k), iMM(f,f);
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
  Eigen::BDCSVD<Eigen::MatrixXf> svd(GEN, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXf V = svd.matrixV();
  Eigen::VectorXf D = svd.singularValues();
  Eigen::MatrixXf X = svd.matrixU() * D.array().matrix().asDiagonal();
  int p = X.cols();
  
  // Sum of squares of X
  if(verb){ Rcpp::Rcout << "Set starting values for coefficients and variances\n"; }
  Eigen::MatrixXf XX(p,k); 
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  Eigen::VectorXf MSX = XX.colwise().sum().array(); MSX=MSX.array()*iN.array();
  // Variances
  iN = (n.array()-f).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSX.array()).matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  Eigen::VectorXf TrDinvXX(k);
  Eigen::MatrixXf Dinv(p,k);
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  // Bending
  float deflate = 1.0, deflateMax = 0.9;
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0;
  // Convergence control
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0;
  int numit = 0;
  // Prior for stability
  float df0 = 1.1;
  Eigen::MatrixXf Sb = vb*df0;
  Eigen::VectorXf Se = ve*df0;
  Eigen::VectorXf iNp = (n.array()+df0-f).inverse();
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
  Eigen::MatrixXf hat = X*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( M * MU.col(i) + hat.col(i)).array(); }
  Eigen::MatrixXf beta = V*b;
  
  // Correlations
  if(verb){ Rcpp::Rcout << "Estimating correlations\n"; }
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  if(verb){  Rcpp::Rcout << "Convergence statistics\n"; }
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
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
  float tmp, VarE, dff, MSB;
  Eigen::MatrixXf ZM(n0,f), PZM(n0,f), ZPZ(f,f);
  Eigen::MatrixXf QtlVar(m,k), QtlHer(m,k), Gamma(m,f), Wald(m,k), PopWald(m,f);
  Eigen::VectorXf mum(f), ZPy(f), gamma(f), iC(n0), tmpCoef(n0);
  Eigen::VectorXf yZg(n0), PyZg(n0);
  Eigen::MatrixXf MSEx(m,k);
  
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
SEXP MLM(Eigen::MatrixXf Y, Eigen::MatrixXf X, Eigen::MatrixXf Z,
         int maxit = 500, float logtol = -8, int cores = 1){
  
  // Basic info
  float df0 = 1.1; 
  if(cores!=1) Eigen::setNbThreads(cores);
  int k = Y.cols(), n0 = Y.rows(), f = X.cols(), p = Z.cols();
  
  // Incidence matrix W
  Eigen::MatrixXf W(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        W(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ W(i,j) = 1.0;}}}
  Eigen::VectorXf n = W.colwise().sum();
  Eigen::VectorXf iN = (n.array()-f).inverse();
  
  // Compute SY, SZ
  Eigen::MatrixXf y(n0,k), WX(n0,f), MU(f,k), iXX(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WX.col(j)=X.col(j).array()*W.col(i).array();}
    iXX = (WX.transpose()*WX).inverse();
    MU.col(i) = iXX * WX.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WX*MU.col(i) ).array()*W.col(i).array(); }
  iXX = (X.transpose()*X).inverse();
  for(int j=0; j<p; j++){ Z.col(j) = (Z.col(j) - X*(iXX*X.transpose()*Z.col(j))).array(); }
  
  // Sum of squares of Z
  Eigen::MatrixXf ZZ(p,k); 
  for(int i=0; i<p; i++){ ZZ.row(i) = Z.col(i).array().square().matrix().transpose() * W;}
  Eigen::VectorXf TrZSZ = ZZ.colwise().sum().array();
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // Variances
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/ (TrZSZ.array()*iN.array())  ).matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXf tilde = Z.transpose() * y;
  
  // Bending
  float deflate = 1.0, deflateMax = 0.9;
  Eigen::MatrixXf A = vb*1.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0;
  
  // Prior for stability
  Eigen::MatrixXf Sb = vb*df0;
  Eigen::VectorXf Se = ve*df0;
  Eigen::VectorXf iNp = (n.array()+df0-f).inverse();
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  Eigen::MatrixXf beta0(p,k);
  float cnv = 10.0;
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
  Eigen::MatrixXf hat = Z*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( X * MU.col(i) + hat.col(i)).array(); }
  
  // Genetic correlations
  Eigen::MatrixXf GC(k,k);
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

// [[Rcpp::export]]
SEXP MRR(Eigen::MatrixXf Y,
         Eigen::MatrixXf X,
         int maxit = 500,
         float tol = 10e-9,
         int cores = 1,
         bool TH = false,
         float NonLinearFactor = 0.0,
         bool InnerGS = false,
         bool NoInversion = false,
         bool HCS = false,
         bool XFA = false,
         int NumXFA = 2,
         float prior_R2 = 0.5,
         float gc_prior_df = 0.5, 
         float var_prior_df = 0.0, 
         float weight_prior_h2 = 0.0,
         float weight_prior_gc = 0.0,
         float PenCor = 0.0,
         float MinCor = 1.0,
         float uncorH2below = 0.0,
         float roundGCupFrom = 1.0,
         float roundGCupTo = 1.0,
         float roundGCdownFrom = 1.0,
         float roundGCdownTo = 0.0,
         float bucketGCfrom = 1.0,
         float bucketGCto = 1.0,
         float DeflateMax = 0.9,
         float DeflateBy = 0.005,
         bool OneVarB = false,
         bool OneVarE = false,
         bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXf xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXf ve = vy * (1-prior_R2);
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  Eigen::VectorXf vbInit = ((vy*prior_R2).array()/MSx.array());
  Eigen::VectorXf veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  float tmp;
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
  Eigen::MatrixXf tilde = X.transpose() * y;
  Eigen::VectorXf TrDinvXSX(k);
  Eigen::MatrixXf Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXf Sb = vb*var_prior_df;
  Eigen::VectorXf Se = ve*var_prior_df;
  Eigen::VectorXf iNp = (n.array()+var_prior_df-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  float bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0;
  int numit = 0;
  float logtol = log10(tol);
  
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
  Eigen::MatrixXf W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXf iVeWj = iVe*1.0;
  Eigen::VectorXf tmpW(p);
  float maxW, minW;
  
  // Objects for other variance structures
  float gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(vb);
  Eigen::MatrixXf UDU(k,k);
  
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
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc_prior_df*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
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
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
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
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Resize convergence vectors
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
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
