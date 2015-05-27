\name{pathway.dat}
\alias{pathway.dat}
\docType{data}
\title{
gene expression data set for pathway selection
}
\description{
p53 cancer cell line gene expression data set
}
\usage{data(pathway.dat)}
\format{
A list of two variables.
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{expression}}{a 50-by-4301 matrix that stores the gene expression data}
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
head(pathway.dat$pathways)
}
\keyword{datasets}
