\name{makeInpDb}

\alias{makeInpDb}

\title{
  Builds an individual DB from online files at inparanoid that is in
  turn meant for generating Inparanoid 8 based objects.
}
\description{
  This is really meant to be used by AnnotationHubData for eventual
  exposure through the AnnotationHub.  Users who are interested in
  Inparanoid 8 should look for the relevant objects to be in
  AnnotationHub.  This function is just a database builder that allows
  us to put the data up there.  So this is not really intended for use
  by end users.
}
\usage{
  makeInpDb(dir, dataDir)
}
\arguments{  
  \item{dir}{the full path to the online Inparanpoid 8 resource.}
  \item{dataDir}{directory where the source inparanoid.sqlite DB is}
}
\value{ A SQLite DB is produced but this is not returned from the
  function. loadDb is expected to be used by AnnotationHub to convert
  these into objects after the sqlite DB is downloaded from the hub server. }
\author{
  M. Carlson 
}
\examples{
## There are paths to resource data for each set of files at Inparanoid.
## This is how you can turn those files inta sqlite DB.
if(interactive()){
 db <- makeInpDb(
        dir="http://inparanoid.sbc.su.se/download/current/Orthologs/A.aegypti/",
        dataDir=tempdir())
}
}
