\name{@ANNOBJPREFIX@WORMBASE}
\alias{@ANNOBJPREFIX@WORMBASE}
\title{Map Wormbase accession numbers with Entrez Gene identifiers}  
\description{ 
  @ANNOBJPREFIX@WORMBASE is an R object that contains mappings between 
  Entrez Gene identifiers and Wormbase accession numbers.    
}

\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} to Wormbase
  Accession Numbers.

  Mappings were based on data provided from ensembl
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
x <- @ANNOBJPREFIX@WORMBASE
# Get the entrez gene IDs that are mapped to a Wormbase ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Wormbase gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}

\keyword{datasets}



