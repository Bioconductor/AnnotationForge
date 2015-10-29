\name{@ANNOBJPREFIX@ZFIN}
\alias{@ANNOBJPREFIX@ZFIN}
\alias{@ANNOBJPREFIX@ZFIN2PROBE}
\title{Map Zfin accession numbers with Entrez Gene identifiers}  
\description{ 
  @ANNOBJPREFIX@ZFIN is an R object that contains mappings between 
  Entrez Gene identifiers and Zfin accession numbers.    
}

\details{
  This object is a simple mapping of manufacturer identifiers to Zfin
  Accessions.

  Mappings were based on data provided by: @ZFINSOURCE@
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
x <- @ANNOBJPREFIX@ZFIN
# Get the entrez gene IDs that are mapped to an Zfin ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Zfin IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}

\keyword{datasets}



