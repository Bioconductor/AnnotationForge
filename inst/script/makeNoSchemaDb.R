## test data lets load up a couple tab files from /extdata...
## file finch_info.txt now has data I can use (for examples)
library(AnnotationForge)

finchFile <- system.file("extdata","finch_info.txt",package="AnnotationForge")
finch <- read.table(finchFile,sep="\t")

## not that this is how it should always be, but that it *could* be this way.
fSym <- finch[,c(2,3,9)]
fSym <- fSym[fSym[,2]!="-",]
fSym <- fSym[fSym[,3]!="-",]
colnames(fSym) <- c("GENEID","SYMBOL","GENENAME")

fChr <- finch[,c(2,7)]
fChr <- fChr[fChr[,2]!="-",]
colnames(fChr) <- c("GENEID","CHROMOSOME")


finchGOFile <- system.file("extdata","GO_finch.txt",package="AnnotationForge")
fGO <- read.table(finchGOFile,sep="\t")
fGO <- fGO[fGO[,2]!="",]
fGO <- fGO[fGO[,3]!="",]
colnames(fGO) <- c("GENEID","GO","EVIDENCE")

## Now make a list
data <- list(gene_info=fSym, chromosome=fChr, go=fGO)
genus <- "Taeniopygia"
species <- "guttata"
dbName <- AnnotationForge:::.generateOrgDbName(genus,species)
## this becomes the file name for the DB
dbfile <- paste(dbName, ".sqlite", sep="")



## Test DB building:
## AnnotationForge:::makeOrgDbFromDataFrames(data, genus, species, dbfile)

# debug(AnnotationForge:::.makeAnnDbPkg)
# debug(AnnotationForge:::.createAnnotPackage)

## or test pkg building
AnnotationForge:::makeOrgPackage(data=data,
                                 version="0.1",
                                 maintainer="Some One <so@someplace.org>",
                                 author="Some One <so@someplace.org>",
                                 outputDir = ".",
                                 tax_id="59729",
                                 genus=genus,
                                 species=species)






## NEXT UP: lets make an actual template in AnnotationDbi so that I
## can start making these things as packages.
