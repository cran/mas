\name{soy}
\alias{y}
\alias{Y}
\alias{Z}
\alias{pop}
\alias{map}
\title{ Soybean dataset }
\description{

Toy dataset part of the SoyNAM. Plant height observed in a few soybean families in multiple environements.

}

\details{

Data was obtained using the \code{BLUP} from the R package SoyNAM. SNPs are coded as 012, where 0 = homozygous as hub parent, 1 = heterozygous, 2 = homozygous as the founder parent. Details of the soybean nested association mappging (SoyNAM) population are provided by Song et al. (2017), Xavier et al. (2018) and Diers et al. (2018).

}
\references{

Song, Q., Yan, L., Quigley, C., Jordan, B.D., Fickus, E., Schroeder, S., Song, B.H., Charles An, Y.Q., Hyten, D., Nelson, R. and Rainey, K., 2017. Genetic characterization of the soybean nested association mapping population. The Plant Genome, 10(2), pp.plantgenome2016-10.

Xavier, A., Jarquin, D., Howard, R., Ramasubramanian, V., Specht, J.E., Graef, G.L., Beavis, W.D., Diers, B.W., Song, Q., Cregan, P.B. and Nelson, R., 2018. Genome-wide analysis of grain yield stability and environmental interactions in a multiparental soybean population. G3: Genes, Genomes, Genetics, 8(2), pp.519-529.

Diers, B.W., Specht, J., Rainey, K.M., Cregan, P., Song, Q., Ramasubramanian, V., Graef, G., Nelson, R., Schapaugh, W., Wang, D. and Shannon, G., 2018. Genetic architecture of soybean yield and agronomic traits. G3: Genes, Genomes, Genetics, 8(10), pp.3367-3375.

}
\examples{
data(soy)
}