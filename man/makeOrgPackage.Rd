\name{makeOrgPackageFromNCBI}

\alias{makeOrgPackageFromNCBI}
% \alias{makeOrgPackageFromBiomart}


\title{
  Making an organism package from annotations available from
  NCBI.
}
\description{
  The \code{makeOrgPackageFromNCBI} function allows the user to make an
  organism package from NCBI annotations available from the NCBI.
}
\usage{

  makeOrgPackageFromNCBI(
    version=,
    maintainer,
    author,
    outputDir=getwd(),
    tax_id,
    genus=NULL,
    species=NULL,
    NCBIFilesDir=getwd(),
    databaseOnly=FALSE,
    useDeprecatedStyle=FALSE,
    rebuildCache=TRUE)
}

\arguments{
  \item{version}{What is the version number for this package?}
  \item{maintainer}{Who is the package maintainer? (must include email to be 
    valid)}
  \item{author}{Who is the creator of this package?}
  \item{outputDir}{A path where the package source should be assembled.}
  \item{tax_id}{The Taxonomy ID that represents your organism.  (NCBI has
  a nice online browser for finding the one you need)}
  \item{genus}{Single string indicating the genus.}
  \item{species}{Single string indicating the species.}
  \item{NCBIFilesDir}{If you want the files to be saved locally that
  were used to make your DB, then specify a directory here.  Otherwise
  leave this NULL.  This argument is mostly useful for troubleshooting
  since many of the files are large.  Most users should ignore it.}
  \item{databaseOnly}{Only generate a DB without a package (TRUE or FALSE)}
  \item{useDeprecatedStyle}{For legacy support in those cases where it
  will actually work.  Should we use the older database and package
  style (with bimaps)? (TRUE or FALSE) Also: this can't work with
  unsupported organisms}
  \item{rebuildCache}{Logical. When TRUE, the cache of local files is 
  refreshed each day. When FALSE, the package is created from existing files 
  in \code{outputDir} and nothing is downloaded. Used internally
  and for testing.}
}

\details{
  \code{makeOrgPackageFromNCBI} downloads multiple files and assembles a 12 GB
  database in \code{NCBIFilesDir}. The first time the function is run it may
  take over an hour; subsequent calls reuse files from the cache and are
  much faster. The default behavior of \code{makeOrgPackageFromNCBI} attempts to 
  refresh the cached files each day (suppress with \code{rebuildCache = FALSE}).
}

\value{Nothing returned to the R session. Just creates an organism
annotation package.}

\author{
  M. Carlson
}


\examples{

\dontrun{
## Makes an organism package for Zebra Finch from NCBI:

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "59729",
                       genus = "Taeniopygia",
                       species = "guttata")

}

}
