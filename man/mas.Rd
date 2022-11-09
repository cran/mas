\name{mas-package}
\alias{mas-package}
\alias{mas}
\docType{package}
\title{\packageTitle{mas}}
\description{\packageDescription{mas}}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{mas}
\packageIndices{mas}
}
\author{
\packageAuthor{mas}
Maintainer: \packageMaintainer{mas}
}
\references{

Wei, J. and Xu, S., 2016. A random-model approach to QTL mapping in multiparent advanced generation intercross (MAGIC) populations. Genetics, 202(2), pp.471-486.

Xavier, A. and Habier, D., 2022. A new approach fits multivariate genomic prediction models efficiently. Genetics Selection Evolution, 54(1), pp.1-15.

}
\examples{

# load the toy dataset
data( soy ) 

# run gwas
fit1 = gwas(y, Z, pop) 
# adjust variances
fit2 = CorrectBeavis( fit1 ) 

# Compare before and after correction
plot( fit1, h2=TRUE, col=8, pch=20) # display QTL h2
plot( fit2, h2=TRUE, add=TRUE, pch=20, type='o') # adjusted QTL h2
legend('topleft',pch=16,col=c(8,1),c('Before correction','After Beavis correction'))


}
