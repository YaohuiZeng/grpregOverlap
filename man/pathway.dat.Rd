\name{pathway.dat}
\alias{pathway.dat}
\docType{data}
\title{
Gene expression and pathway information of p53 cancer cell lines
}

\description{
The data set contains gene expression data and pathway (group) information of p53 
cancer cell lines. The mutational status of the p53 gene for 50 cell lines is 
recorded, with 17 classified as normal and 33 as carrying mutations. Pathway 
information of the genes are from the C2 catalog of the Initial Catalog of Human
Gene Sets, or MSigDB 1.0 (Subramanian et al., 2005).
}

\usage{data(pathway.dat)}
\format{
The raw data files of gene expression and pathway information can be found via
links in \strong{Source} section below. The raw data is preprocessed such that
only 308 pathways with size between 15 and 500 are included. Then 4301 genes in 
those selected pathways are chosen.

A list of three variables included in \code{pathway.dat}:
  \describe{
    \item{\code{expression}}{
    a 50-by-4301 matrix that records the gene expression data. Used as design 
    matrix.
    }
    \item{\code{mutation}}{
    a 1-by-50 binary response vector recording the mutational status: 1 = normal;
    0 = mutation. Used as response vector.
    }
    \item{\code{pathways}}{
    a list of 308 vectors. Each contains the names of genes in that pathway. Used
    as group information.
    }
  }
}

\source{
The raw data files can be downloaded via \url{http://www.broadinstitute.org/gsea/datasets.jsp}.
}

\references{
  \itemize{
    \item Subramanian, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. \emph{Proceedings of the National Academy of Sciences of the United States of America}, \strong{102}(43), 15545-15550. \url{http://www.pnas.org/content/102/43/15545.short}
  }
}

\examples{
data(pathway.dat)
pathway.dat$expression[1:10, 1:10]
pathway.dat$mutation
head(pathway.dat$pathways)
}
\keyword{datasets}
