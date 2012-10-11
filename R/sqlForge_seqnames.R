
###############################################################################
###############################################################################
###############################################################################
###                   Code for generating seqnames.db.                      ### 
###############################################################################
###############################################################################
###############################################################################

###############################################################################
## organisms to start with:
## arabidopsis
## celegans
## drosophila
## human
## mouse
## rat
## yeast


###############################################################################
## To generate a fresh DB simply:
## 1) place any new csv files into the extdata/dataFiles dir and then:
##    R CMD INSTALL AnnotationDbi
## 2) Call generateSeqnames.db()
##



###############################################################################
## function to help create and populate a seqnames DB from csv files.
## 'specString' is species name (must match a file), and
## 'where' is path to location for the DB

writeTable <- function(specString, con){
  tables <- dbListTables(con)
  if(any(specString %in% tables)){
    stop(paste("there is already a",specString,"table"))
  }
  df <- system.file("seqnames-template","inst","extdata","dataFiles",
                     package="AnnotationForge")
  data <- read.table(file.path(df, paste(specString,".csv",sep="")),
                     header=TRUE, sep="\t",
                     colClasses="character", check.names=FALSE)
  message("Creating table: ",specString)
  sqliteWriteTable(con, specString, value = data, row.names = FALSE)
}



###############################################################################
## regenDB() is for 1) creating the tables, 2) reading in the appropriate
## files, and populating the DB.
## 'where' defines the path to the created DB.

regenDb <- function(where){
  require(RSQLite)
  con <- dbConnect("SQLite", dbname=file.path(where,"seqnames.sqlite"))
  ## list files in the extData/tabfiles dir
  dfs <- system.file("seqnames-template","inst","extdata","dataFiles",
                     package="AnnotationForge")
  files <- gsub(".csv","",dir(dfs))
  ## Then for each file, call makeTable
  lapply(files, writeTable, con=con)
}

## usage example:
## regenDb(where=".")




## A function to create the a new DB from the template.
generateSeqnames.db <- function(version, outdir="."){
  require(RSQLite)
  ## use creatPackage (from Biobase)
  symbolValues <- list(VERSION=version,
                       ANNOTATIONDBIVERSION = "1.19.4")
  message("Creating the package directories")
  createPackage(pkgname="seqnames.db",
                destinationDir = outdir,
                originDir = system.file("seqnames-template",
                  package="AnnotationForge"),
                symbolValues = symbolValues,
                unlink = TRUE,
                quiet = TRUE)
  
  message("Creating the database")
  ## now generate the DBs
  regenDb(where=file.path(getwd(),"seqnames.db","inst","extdata"))
}


## require(AnnotationForge); AnnotationForge:::generateSeqnames.db(version="1.1.0")


## TODO: export generateSeqnames.db and write a manual page for it.



## ## TODO fix this bug!:
##   library("TxDb.Hsapiens.UCSC.hg19.knownGene")
##   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
##   species = species(txdb)
##   style = "NCBI"
##   listAllSupportedStylesBySpecies(species)

