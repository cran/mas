
CorrectBeavis = function(fit,verb=TRUE){
  
  # Check class
  if(!is(fit,"GenModel")) stop("Input object must be of class GenModel.")
  
  # Setting threshold
  m = nrow(fit$GWAS$Wald)
  alpha = 0.05/m
  q = ncol(fit$GWAS$Membership)
  tw = qchisq(1-alpha,df=q)
  if(verb) cat('Bonferroni-adjusted Wald Threshold (tw) =',round(tw,2),'\n')
  
  # Shizhong's function
  correct = function(w,mse){
    a<-tw
    b<-Inf
    beta<-(tw+2)/(q+tw+2)
    fn<-function(d,w) ifelse(d<100,extrunc("chisq",a=tw,ncp=d,df=q)-w,d+q-w)
    ncp_tm<-function(w){
      f<-extrunc("chisq",a=tw,df=q,ncp=0)
      if(f>w){ delta = 0 }else{
        delta<-uniroot(f=fn,w=w,lower=0,upper=1000,tol=1e-5)$root}
      return(delta)}
    xbeta<-function(x) ncp_tm(x)-beta*(x-tw)
    ncp_m<-function(w){
      J0<-extrunc("chisq",a=tw,df=q,ncp=0)
      xb<-min(uniroot(f=xbeta,lower=J0,upper=100,tol=1e-5)$root)
      if(a<=w && w<xb){ delta<-beta*(w-a)}else if(w<a){
        delta<-0}else{delta<-ncp_tm(w)}
      return(delta)}
    if(w>tw){ delta<-ncp_m(w=w); vqtl.correct<-delta*mse/(n-q) }else{ vqtl.correct<-0 }   
    return(vqtl.correct)
  }
  
  # Adjustment
  if(verb) cat('Correcting Var(QTL) and H2(QTL)\n')
  k = ncol(fit$GWAS$Wald)
  CorrectedVarQTL = matrix(NA,m,k)
  for(i in 1:k){
    cat('Trait',i,'\n')
    mse=fit$GWAS$MSE[,i]
    wald=fit$GWAS$Wald[,i]
    n=fit$GWAS$N[i]
    CorrectedVarQTL[,i] = mapply(correct,w=wald,mse=mse)
  }
  
  # Updating object
  if(verb) cat('All done!\n')
  fit2 = fit
  CorrectedVarQTL[CorrectedVarQTL<0] = 0
  fit2$GWAS$QTLvar = CorrectedVarQTL
  qtlh2 = CorrectedVarQTL/(fit$GWAS$MSE+CorrectedVarQTL)
  qtlh2[qtlh2<0] = 0
  fit2$GWAS$QTLh2 = qtlh2
  fit2$GWAS$BeavisCorrection = TRUE
  
  # Output
  return(fit2)
  
}
