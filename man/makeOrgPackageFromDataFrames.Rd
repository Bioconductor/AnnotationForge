\name{makeOrgPackage}

\alias{makeOrgPackage}


\title{
  Making an organism package from annotations available from a set of named
  data.frames.
}
\description{
  The \code{makeOrgPackage} function allows the user to make an
  organism package from any collection of data frames that are united by
  a common gene ID.
}
\usage{

  makeOrgPackage(...,
                 version,
                 maintainer,
                 author,
                 outputDir = ".",
                 tax_id,
                 genus,
                 species,
                 goTable=NA)
}

\arguments{
  \item{...}{A set of data.frames containing annotation data.  Each of
  these arguments must be named.  Those names will become the names of
  the tables in the final database.  Also, there are no rownames for
  these data.frames, and the colnames are the names that will be used as
  extractable fields in the final package.  In other words they will be
  what comes back when you call \code{columns} and \code{keytypes}
  Finally, the 1st column of every data.frame must be labeled GID, and
  correspond to a gene ID that is universal for the entire set of
  data.frames.  The GID is how the different tables will be joined
  internally}
  \item{version}{What is the version number for this package? format: \'x.y.z\'}
  \item{maintainer}{Who is the package maintainer? (must include email to be 
    valid)}
  \item{author}{Who is the creator of this package?}
  \item{outputDir}{A path where the package source should be assembled.}
  \item{tax_id}{The Taxonomy ID that represents your organism.  (NCBI has
  a nice online browser for finding the one you need)}
  \item{genus}{Single string indicating the genus.}
  \item{species}{Single string indicating the species.}
  \item{goTable}{By default, this is NA, but if one of your '...'
  data.frames has GO annotations, then this name will be the name of
  that argument.  When you specify this, \code{makeOrgPackage} will
  process that data.frame to remove extra GO terms (that are too new for
  the current GO.db) and also will generate a table for GOALL data
  (based on ancestor terms for each mapping from GO.db).  This table
  will also be checked to make sure that it has exactly THREE columns,
  that must be named GID, GO and EVIDENCE.  These must correspond to the
  gene IDs, GO IDs and evidence codes respectively.  GO IDs should be
  formatted like this to work with other DBs in the project:
  \'GO:XXXXXXX\'. }
}


\value{The path to the package that just created. This is useful for
  calling install.packages as the next step.}

\author{
  M. Carlson
}


\examples{

if(interactive()){

## Makes an organism package for Zebra Finch data.frames:
finchFile <- system.file("extdata","finch_info.txt",package="AnnotationForge")
finch <- read.table(finchFile,sep="\t")

## not that this is how it should always be, but that it *could* be this way.
fSym <- finch[,c(2,3,9)]
fSym <- fSym[fSym[,2]!="-",]
fSym <- fSym[fSym[,3]!="-",]
colnames(fSym) <- c("GID","SYMBOL","GENENAME")

fChr <- finch[,c(2,7)]
fChr <- fChr[fChr[,2]!="-",]
colnames(fChr) <- c("GID","CHROMOSOME")

finchGOFile <- system.file("extdata","GO_finch.txt",package="AnnotationForge")
fGO <- read.table(finchGOFile,sep="\t")
fGO <- fGO[fGO[,2]!="",]
fGO <- fGO[fGO[,3]!="",]
colnames(fGO) <- c("GID","GO","EVIDENCE")

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="59729",
               genus="Taeniopygia",
               species="guttata",
               goTable="go")

## then you can call install.packages based on the return value
install.packages("./org.Tguttata.eg.db", repos=NULL)

}

}
