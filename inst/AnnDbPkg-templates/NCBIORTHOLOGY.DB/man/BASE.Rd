\name{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@.db}
\alias{@ANNOBJPREFIX@}
\title{Bioconductor annotation data package}
\description{

  Welcome to the @ANNOBJPREFIX@.db annotation Package. The purpose of
  this package is to provide orthology mappings between species, based
  on NCBI Gene IDs and NCBI orthology mappings. This package is updated biannually.

  Objects in this package are accessed using the \code{select()} interface.
  See ?select in the AnnotationDbi package for details.

  Please note that this package is slightly different from other
  annotation packages; the keytypes and columns are species names,
  separated by a period (e.g. Homo.sapiens), and the keys are the NCBI
  Gene IDs that are to be mapped between species. See examples for more details.
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
## Map NCBI Gene IDs from human to mouse
select(Orthology.eg.db, as.character(c(1,9,13:22)),"Mus.musculus","Homo.sapiens")

## get list of available species
head(keytypes(Orthology.eg.db))

## get list of Gene IDs for Electric eel
eeKeys <- keys(Orthology.eg.db, "Electrophorus.electricus")
head(eeKeys)

## map all availble Electric eel genes to human
ee2human <- select(Orthology.eg.db, eeKeys, "Homo.sapiens", "Electrophorus.electricus")

}
\keyword{datasets}

