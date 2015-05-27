\name{pathway.dat}
\alias{pathway.dat}
\docType{data}
\title{
p53 cancer cell line gene expression data
}
\description{
gene expression data set for pathway selection
}
\usage{data(pathway.dat)}
\format{
A list of three variables.
  \describe{
    \item{\code{expression}}{a 50-by-4301 matrix that stores the gene expression data}
    \item{\code{mutation}}{a 1-by-50 binary response vector}
    \item{\code{pathways}}{a list of 308 vectors. Each contains the names of genes that in the pathway}
  }
}
\details{
See GSEA paper for details. Only include pathways with size between 15 and 500, and genes in those pathways.
}
\source{
GSEA: \url{http://www.pnas.org/content/102/43/15545.full?ref=Guzels.TV}
}
\references{
GSEA: \url{http://www.pnas.org/content/102/43/15545.full?ref=Guzels.TV}
}

\examples{
data(pathway.dat)
pathway.dat$expression[1:10, 1:10]
pathway.dat$mutation
head(pathway.dat$pathways)
}
\keyword{datasets}
