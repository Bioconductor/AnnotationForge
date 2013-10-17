## test to script how I expect to call makeChipPackage
library(AnnotationForge)

finchFile <- system.file("extdata","finch_info.txt",package="AnnotationForge")
finch <- read.table(finchFile,sep="\t")
geneIds <- finch$V2
probeNames <- paste("probe", 1:length(geneIds), sep="")
probeFrame <- data.frame(probes=probeNames, genes=geneIds)

orgPkgName <- "org.TguttataTestingSubset.eg.db"
genus <- "Taeniopygia"
species <- "guttataTestingSubset"
dbName <- AnnotationForge:::.generateOrgDbName(genus,species)
## this becomes the file name for the DB
dbfile <- paste(dbName, ".sqlite", sep="")
tax_id <- "59729"

## Test DB building:
AnnotationForge:::makeOrgDbFromDataFrames(probeFrame, orgPkgName, tax_id,
                                          genus, species, dbfile)

# debug(AnnotationForge:::.makeAnnDbPkg)
# debug(AnnotationForge:::.createAnnotPackage)


## or test pkg building
AnnotationForge:::makeOrgPackage(probeFrame,
                                 orgPkgName,
                                 version="0.1",
                                 maintainer="Some One <so@someplace.org>",
                                 author="Some One <so@someplace.org>",
                                 outputDir = ".",
                                 tax_id="59729",
                                 genus=genus,
                                 species=species)

