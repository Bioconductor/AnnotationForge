\name{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@}
\title{Bioconductor annotation data package}
\description{
  
  Welcome to the @ANNOBJPREFIX@.db annotation Package.  The purpose of
  this package is to provide detailed information about the
  @ANNOBJPREFIX@ platform. Bioconductor attempts to update this package
  biannually with the most recently available inparanoid data, but
  inparanoid data sources have not historically updated that frequently.
  Please check inparanoid to learn more about when they have last
  updated their databases.
  
  Objects in this package are accessed using the \code{select()} interface.
  See ?select in the AnnotationDbi package for details.
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          \code{keys()}, \code{columns()} and \code{select()}.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.
columns(@ANNOBJPREFIX@.db)

## Bimap interface:
## The 'old style' of interacting with these objects is manipulation as 
## bimaps. While this approach is still available we strongly encourage the
## use of select().
ls("package:@ANNOBJPREFIX@.db")
}
\keyword{datasets}

