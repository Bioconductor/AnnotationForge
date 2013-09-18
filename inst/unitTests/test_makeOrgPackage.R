## tests for helpers for makeOrgPackage

## 1st prep some data for testing.
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

## Now make a list (for testing internal stuff)
data <- list(gene_info=fSym, chromosome=fChr, go=fGO)
## And set up come other variables
genus <- "Taeniopygia"
species <- "guttata"
goTable <- "go"
maintainer <- "Some One <so@someplace.org>"
author <- "Some One <so@someplace.org>"
outputDir <- tempdir()
dbName <- AnnotationForge:::.generateOrgDbName(genus,species)
dbFileName <- file.path(outputDir,paste0(dbName, ".sqlite"))
tax_id <- "59729"
goTable <- "go"

## tests

## Mostly this will either break or it will work, but there are a few
## things that should be checked...

test_.generateOrgDbName <- function(){
    dbName <- AnnotationForge:::.generateOrgDbName(genus,species)
    checkTrue(dbName=="org.Tguttata.eg")
}

test_addOntologyData <- function(){
    goData <- data[[goTable]]
    res <- AnnotationForge:::.addOntologyData(goData)
    checkTrue(dim(res)[1] == dim(goData)[1])
    checkTrue(dim(goData)[2] == 3)
    checkTrue(all(names(res) == c('GID','GO','EVIDENCE','ONTOLOGY')))
}

test_expandGOFrame <- function(){
    goData <- data[[goTable]]
    goData <- goData[goData[["GO"]] %in% Lkeys(GOTERM),]
    goData <- AnnotationForge:::.addOntologyData(goData)
    gbp <- goData[goData$ONTOLOGY=="BP",c("GID","GO","EVIDENCE")]
    names(gbp) <- c("gene_id","go_id","evidence")
    require(GO.db)
    res <- AnnotationForge:::.expandGOFrame(gbp, GOBPANCESTOR)
    checkTrue(dim(res)[1] > dim(goData)[1]) ## expansion must have happened
    checkTrue(all(names(res) == c('gene_id','go_id','evidence')))
}


test_makeOrgPackage <- function(){
    pkgNm <- AnnotationForge:::makeOrgPackage(gene_info=fSym,
                                       chromosome=fChr,
                                       go=fGO,
                                       version="0.1",
                                       maintainer="Some One <so@someplace.org>",
                                       author="Some One <so@someplace.org>",
                                       outputDir = outputDir,
                                       tax_id="59729",
                                       genus=genus,
                                       species=species,
                                       goTable="go")

    
    checkTrue(pkgNm==file.path(outputDir,"org.Tguttata.eg.db"))
##     ## Then install it.
##     install.packages(pkgNm, repos=NULL)

##     ## Now test the output of select, cols, keytypes, and keys...
##     library(org.Tguttata.eg.db)
##     ## debug(AnnotationDbi:::.noSchemaCols)
##     columns(org.Tguttata.eg.db)
##     ## so that change will also work for keytypes
##     keytypes(org.Tguttata.eg.db)

    
    ## finish by removing it
##     remove.packages("org.Tguttata.eg.db")

    ## the install and remove options don't work with R CMD check :(
    ## ALSO the keytypes() etc.  calls should be checked in AnnotationDbi
}

