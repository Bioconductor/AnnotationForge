\name{@ANNOBJPREFIX@GENETYPE}
\alias{@ANNOBJPREFIX@GENETYPE}
\title{Map between Entrez Gene Identifiers and Gene Type}
\description{
  @ANNOBJPREFIX@GENETYPE is an R object that provides mappings between
  entrez gene identifiers and gene types. 
}
\details{
  Each entrez gene identifier is mapped to the reported type of 
  the corresponding gene. An \code{NA} is reported if there is no known
  type for a given gene.
  
}

\references{
  \url{https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene}
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
x <- @ANNOBJPREFIX@GENETYPE
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
# Get the entrez gene identifiers that are mapped to a gene type
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the entrez gene ID for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}

\keyword{datasets}


