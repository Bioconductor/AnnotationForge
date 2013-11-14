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




## TODO: you forgot to add accesions tables.  You really just only
## need this for the legacy chip format.  It's a two column table that
## looks like this:
## CREATE TABLE accessions (probe_id VARCHAR(80),accession VARCHAR(20));
## CREATE INDEX Fgbprobes ON accessions (probe_id);
## And it should be created whenever .makeLegacyProbesTable() is called.
## I don't want that to be part of the new schema-less ChipDbs though.
## Checking: For the main function, it can just be an option that you
## can pass in if you are using a legacy package.  (A warning will
## have to be issued if you use it with a schema-less org package).


## TODO: test to make sure that I can do this with a probeFrame where
## some probes have no gene IDs. = OK

## BUT: There is still an oddity, it seems that I need to put entrez
## gene IDs into accessions table too?  That's really strange for
## users...  Maybe the argument should just be left out?



##########################################################################
## Notes on the planned behavior of makeChipPackage()
## the internal table name should be called probes with fields
## (PROBE and GID) to go with new NOSCHEMA style of database.

## What if someone wants to use an old style org package with a
## new style of chip package?  - That situation could have gotten
## complicated BUT chip packages DO specify the org package that
## they are supposed to depend on.

## Best idea I think: based on the org package, I need to generate
## either a NOSCHEMACHIP OR a CHIPDB (there is very little to do
## in either case), so that the user can specify what they want
## and get a package built for their needs.
## SO by default get them a NOSCHEMACHIP.DB, and otherwise it's
## one of the following: NCBICHIP.DB, YEASTCHIP.DB,
## ARABIDOPSISCHIP.DB.  If its one of those "other three" then
## dispatch will be handled through .legacySelect()

## TODO: I am going to have to write a DBSCHEMA detection and dispatch
## helper for finding if it's ARABIDOPSIS, YEAST, NOSCHEMA OR
## something else.
