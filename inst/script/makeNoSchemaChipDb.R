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


## specifics for the DB test
prefix = "fakeChip"
outputDir="."
dbFileName <- file.path(outputDir,paste0(prefix, ".sqlite"))
tax_id <- "59729"

## Test DB building:
AnnotationForge:::makeChipDbFromDataFrame(probeFrame, orgPkgName, tax_id,
                                          genus, species, dbFileName)

# debug(AnnotationForge:::.makeAnnDbPkg)
# debug(AnnotationForge:::.createAnnotPackage)


## or test pkg building
AnnotationForge:::makeChipPackage(prefix=prefix,
                                  probeFrame=probeFrame,
                                  orgPkgName=orgPkgName,
                                  version="0.1",
                                  maintainer="Some One <so@someplace.org>",
                                  author="Some One <so@someplace.org>",
                                  outputDir=outputDir,
                                  tax_id=tax_id,
                                  genus=genus,
                                  species=species)


## next up:

## -make NOSCHEMACHIP.db - done
## -deal with creating multiple "kinds" (old and new) - do we need a
## new kind? - see notes - we do want that, but we have to also have
## legacy support... - done?
## - have to add support for legacy chip mappings (four were never defined - so try to make that work)


## -make changes to AnnotationDbi to support these different package types.

## - have to add support for select methods
