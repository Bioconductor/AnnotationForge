\name{@ANNOBJPREFIX@REJECTORF}
\alias{@ANNOBJPREFIX@REJECTORF}
\title{Rejected Yeast Genes (ORF)}
\description{
  This is based upon Real and rejected yeast ORFs from Kellis et
  al. (2003) to list the rejected genes by this criteria.
}
\details{
  A character vector which contains yeast ORFs which are rejected
  in the reading frame conseration (RFC) test in Kellis et al. (2003).
}

\references{
  Manolis Kellis, Nick Patterson, Matthew Endrizzi, Bruce Birren and
  Eric S. Lander, Sequencing and comparison of yeast species to identify
  genes and regulatory elements.  Nature 423, 241-254 (15 May 2003)
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
head(@ANNOBJPREFIX@REJECTORF)
}

\keyword{datasets}

