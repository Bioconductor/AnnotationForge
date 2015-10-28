\name{@ANNOBJPREFIX@PROTEINGI}
\alias{@ANNOBJPREFIX@PROTEINGI}
\title{Map Entrez Gene IDs to protein GenInfo (GI) numbers}
\description{
 @ANNOBJPREFIX@PROTEINGI is an R object that provides mappings between entrez gene
 identifiers and the protein GI number associated to the gene of interest. 
}
\details{
  Each entrez gene identifier maps to a protein gi number. 
   
  Mappings were based on data provided by: @PROTEINGISOURCE@
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
x <- @ANNOBJPREFIX@PROTEINGI
# Get the entrez gene identifiers that are mapped to a protein GI number
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the PROTEINGI for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}

\keyword{datasets}
