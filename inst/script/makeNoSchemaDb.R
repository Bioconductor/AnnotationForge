## test data lets load up a couple tab files from /extdata...
## file finch_info.txt now has data I can use (for examples)
library(AnnotationForge)

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


## then you can install on the return value
install.packages("./org.Tguttata.eg.db", repos=NULL)


## NEXT UP: lets make an actual template in AnnotationDbi so that I
## can start making these things as packages.
library(org.Tguttata.eg.db)
#debug(AnnotationDbi:::.noSchemaCols)
columns(org.Tguttata.eg.db)
## so that change will also work for keytypes
keytypes(org.Tguttata.eg.db)

## Now for keys I need to read carefully.
## I want to do the same thing for .keys that I did for
library(org.Tguttata.eg.db)
# debug(AnnotationDbi:::.noSchemaKeys)
keys(org.Tguttata.eg.db, "CHROMOSOME")

head(keys(org.Tguttata.eg.db, "GID"))

head(keys(org.Tguttata.eg.db, "SYMBOL", pattern="BDNF"))

head(keys(org.Tguttata.eg.db, "GID", pattern="BDNF", column="SYMBOL"))

head(keys(org.Tguttata.eg.db, "SYMBOL", column="GID"))



## Now I just need select() to work...
library(org.Tguttata.eg.db)
## debug(AnnotationDbi:::.noSchemaSelect)

## TODO: add check to the keys argument to make sure its a character()
select(org.Tguttata.eg.db, keys="100008579", columns="SYMBOL", keytype="GID")
## 
select(org.Tguttata.eg.db, keys="100008579", columns=c("SYMBOL","CHROMOSOME"), keytype="GID")
## This one should warn about the number or rows - it does
select(org.Tguttata.eg.db, keys="100008579", columns="GO", keytype="GID")
## now fixed
select(org.Tguttata.eg.db, keys="100008579", columns=c("GO","EVIDENCE"), keytype="GID")
## What if there is only one table to visit?
select(org.Tguttata.eg.db, keys="BDNF", columns="GENENAME", keytype="SYMBOL")

