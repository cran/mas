gwas = function(y,GEN,m=NULL,...){
  
  if(is.vector(y)){
    Y = matrix(y,ncol=1)
  }else{
    Y=y
  } 
  
  if(is.null(m)){
    M = matrix(rep(1,nrow(Y)))
  }else if(is.vector(m)){
    m = droplevels(factor(m))
    M = model.matrix(~m-1)
  }else if(is.matrix(m)){
    M = m
  } 
  
  fit = GWAS(Y,GEN,M,...)
  return(fit)
}

plot.GenModel = function(x,h2=FALSE,trait=1,add=FALSE,...){
  k = ncol(x$GWAS$Wald)
  m = nrow(x$GWAS$Wald)
  alpha = 0.05/m
  q = ncol(x$GWAS$Membership)
  tw = qchisq(1-alpha,df=q)
  kk = ceiling(sqrt(k))
  if(h2){
    if(add){
      lines(x$GWAS$QTLh2[,trait],...)
    }else{
      plot(x$GWAS$QTLh2[,trait],ylab='QTL heritability',xlab='Marker',...)
    }
  }else{
    if(add){
      lines(x$GWAS$Wald[,trait],...)
    }else{
      plot(x$GWAS$Wald[,trait],ylab='Wald',xlab='Marker',...)
      abline(h=tw,lty=3)
    }
  }
}

print.GenModel = function(x, ...){
  k = ncol(x$GWAS$Wald)
  m = nrow(x$GWAS$Wald)
  alpha = 0.05/m
  q = ncol(x$GWAS$Membership)
  n = round(mean(x$GWAS$N))
  tw = qchisq(1-alpha,df=q)
  cat('Association analysis run on',k,'trait(s),',m,'markers,',n,'individuals from',q,'populations\n')
  h2 = paste(round(x$GBLUP$Heritability,2),collapse = ', ')
  cat('Trait heritability:',h2,' \n')
  as = paste(colSums(x$GWAS$Wald>tw),collapse = ', ')
  cat('Number of significant associations:',as,' \n')
  if("BeavisCorrection"%in%ls(x$GWAS)){cat('Variances were adjusted for the Beavis effect\n')} 
}
