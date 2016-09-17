\name{@ANNOBJPREFIX@ENTREZID}
\alias{@ANNOBJPREFIX@ENTREZID}
\title{Map Systematic ORF identifiers with Entrez Gene identifiers}  
\description{ 
  @ANNOBJPREFIX@ENTREZID is an R object that contains mappings between 
  Systematic ORF accession numbers and NCBI Entrez Gene identifiers.    
}
\details{
  This object is a simple mapping of Entrez Gene identifiers
  \url{https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene} to
  systematic ORF locus Accession Numbers.

  Mappings were based on data provided by: @ENTREZIDSOURCE@
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
x <- @ANNOBJPREFIX@ENTREZID
# Get the ORF IDs that are mapped to an Entrez Gene ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Entrez gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}
\keyword{datasets}



