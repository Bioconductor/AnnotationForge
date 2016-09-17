\name{@ANNOBJPREFIX@ACCNUM}
\alias{@ANNOBJPREFIX@ACCNUM}
\alias{@ANNOBJPREFIX@ACCNUM2EG}
\title{Map Entrez Gene identifiers to GenBank Accession Numbers}  
\description{ 
  @ANNOBJPREFIX@ACCNUM is an R object that contains mappings between 
  Entrez Gene identifiers and GenBank accession numbers.    
}
\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{https://www.ncib.nlm.nih.gov/entrez/query.fcgi?db=gene} to all possible GenBank 
  accession numbers.

  Mappings were based on data provided by: @ACCNUMSOURCE@
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
x <- @ANNOBJPREFIX@ACCNUM
# Get the entrez gene identifiers that are mapped to an ACCNUM
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the ACCNUM for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}

\keyword{datasets}



