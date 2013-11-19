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
## AnnotationForge:::makeChipDbFromDataFrame(probeFrame, orgPkgName, tax_id,
##                                           genus, species, dbFileName)

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


## TESTING select() (for example)
library(fakeChip.db)
keytypes(fakeChip.db)                         ## now works!
columns(fakeChip.db)                          ## now works!
k = head(keys(fakeChip.db))                   ## now works!
k
select(fakeChip.db , k , "SYMBOL","PROBEID")  

## TODO: 1) ORGPKG needs to be in the metadata and then get made by
## template into relevant object. 2) ORGPKG ALSO needs to be inserted
## into the DESCRIPTION file template (as a dependency) - so add a new
## "tag" to the templates.
## BUG:
library(org.TguttataTestingSubset.eg.db)
k = head(keys(org.TguttataTestingSubset.eg.db))
select(org.TguttataTestingSubset.eg.db, k , c("SYMBOL","CHROMOSOME"),"GID")



## test a legacy platform (do a human example)
library(org.Hs.eg.db)
library(AnnotationForge)
geneIds = head(keys(org.Hs.eg.db),n=200)
probeNames <- paste("probe", 1:length(geneIds), sep="")
probeFrame <- data.frame(probes=probeNames, genes=geneIds)

#debug(AnnotationForge:::makeChipPackage)
# debug(AnnotationForge:::.cloneMapMetadata)

AnnotationForge:::makeChipPackage(prefix="fakeHumanChip",
                                  probeFrame=probeFrame,
                                  orgPkgName="org.Hs.eg.db",
                                  version="0.1",
                                  maintainer="Some One <so@someplace.org>",
                                  author="Some One <so@someplace.org>",
                                  outputDir=".",
                                  tax_id="9606",
                                  genus="Homo",
                                  species="sapiens")



## optionally, there should be values for accessions that can be passed in
accs = head(keys(org.Hs.eg.db, keytpye="ACCNUM"),n=300)
probeNames <- paste("probe", 1:length(accs), sep="")
accFrame <- data.frame(probes=probeNames, accs=accs)
## debug(AnnotationForge:::.makeLegacyProbesTable)
AnnotationForge:::makeChipPackage(prefix="fakeHumanChip2",
                                  probeFrame=probeFrame,
                                  orgPkgName="org.Hs.eg.db",
                                  version="0.1",
                                  maintainer="Some One <so@someplace.org>",
                                  author="Some One <so@someplace.org>",
                                  outputDir=".",
                                  tax_id="9606",
                                  genus="Homo",
                                  species="sapiens",
                                  optionalAccessionsFrame=accFrame)


## next up:

## -make NOSCHEMACHIP.db - done
## -deal with creating multiple "kinds" (old and new) - do we need a
## new kind? - see notes - we do want that, but we have to also have
## legacy support... - done?
## - have to add support for legacy chip mappings (four were never defined - so try to make that work)


## -make changes to AnnotationDbi to support these different package types.

## - have to add support for select methods


## We ALMOST can use the 2nd type of package.  (the 1st type will need a bit more work)

## TODO: make sure we support ORGANISM_DB schema for our easy chip packages...

## This now (mostly) works (at least for human legacy pkgs).  The one thing I see wrong is that there are not any ACCNUM values (which makes sense).  So I might have to take steps to stamp out the ACCNUM man page etc. from this kind of package (after the fact)


## for example
keytypes(fakeHumanChip.db)
columns(fakeHumanChip.db)
k = head(keys(fakeHumanChip.db))
select(fakeHumanChip.db , k , "SYMBOL","PROBEID")






##########################################################################
## We have two kinds of chip packages we can make, 1) NOCHIPSCHEMA_DB
## packages ddon't have the same structure internally as classic chip
## packages, but they do support select. 2) classic ones point to a
## known (supported) org package and also will have bimaps and most of
## what is in a classic chip package (map counts are gone though).
## BUT overall, I think these are still BETTER for me because you
## don't need a .db0 to make one.
## There is an extra argument for the 2nd case, for the case where
## users want to have a functioning ACCNUM bimap (and have ACCNUMS).
## But it's optional.  This second kind of package really exists so
## that I can make classic chip packages in this new way (still only
## for supported packages).  If I go this route, I should probably add
## some code to fill in the map_counts.


## TODO: Add better argument checking to make sure columns and
## keytypes are always legitimate! - DONE.
library(org.Hs.eg.db)
k = head(keys(org.Hs.eg.db))
select(org.Hs.eg.db , k , "SYMBOL","ENTREZID")









## TODO:
## This ALSO means that PROBEID is now officially a reserved word for
## making NOSCHEMA org packages.  Also forbidden is the use of probes
## as a table name
## ALSO: genes should be reserved as a table name (and should be already)



