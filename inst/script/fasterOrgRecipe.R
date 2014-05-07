## The aim here is to sketch out and then implement a strategy for
## making all the org packages quickly by calling the makeOrgPackage()
## on code processed from the FTP file downloaded by the NCBI_ftp.R
## source.

## GOALS:
## 1) I want this to be performant.  It cannot take 30 minutes to
## download the whole FTP file and split that every time (and most of
## the time is in the splitting).  I need to do the initial processing
## step once and then run the recipe on some generated stuff that
## comes out of that. - DONE

## 2) I want to pre-filter based on a. whether the tax ID is legit,
## b. whether there are records for that id and c. whether there are
## GO annotations for that organism either at NCBI OR at blast2GO.

## 3) Some examples (see below) don't work, if it's because there is
## no data- FINE.  But the code needs to be more robust AND I need to
## not have these critters in the pre-filtered set.

## 4) It might makes sense to use the makeOrgPackage code as a base
## for making the classic org packages FIRST (it's basically a good
## opportunity to refactor).



## So lets explore an example and see what is actually the slow part...
## slow step is the part where we discard data.


debug(makeOrgDbFromNCBI)
debug(AnnotationForge:::.makeBaseDBFromDLs)
debug(AnnotationForge:::.saveFiles)
debug(AnnotationForge:::.downloadData)
debug(AnnotationForge:::.getFiles)
debug(AnnotationForge:::.writeToNCBIDB)
debug(AnnotationForge:::.indexTaxIds)


## test for older stuff
library(AnnotationForge)
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "9606",
                       genus = "Homo",
                       species = "sapiens",
                       NCBIFilesDir=".",
                       useDeprecatedStyle=TRUE)
## WORKS
debug(AnnotationForge:::prepareDataFromNCBI)

## Another test for building the new data.frames out


debug(AnnotationForge:::prepareDataFromNCBI)


debug(AnnotationForge:::.getBlast2GO)


library(AnnotationForge)

## So this is an example where there are not any refseq or genbank IDs
## (tax ID is too general)
makeOrgPackageFromNCBI(version="0.11",
                       maintainer="g <g.ye@irri.org>",
                       author="g <g.ye@irri.org>",
                       outputDir=".",
                       tax_id="4530",
                       genus="OryZa",
                       species="sativa",
                       NCBIFilesDir=".")

## Whereas this one works just fine
makeOrgPackageFromNCBI(version="0.11",
                       maintainer="g <g.ye@irri.org>",
                       author="g <g.ye@irri.org>",
                       outputDir=".",
                       tax_id="39947",
                       genus="OryZa",
                       species="sativa.japonica",
                       NCBIFilesDir=".",
                       databaseOnly=TRUE)




## The following all works:
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "192222",
#                       genus = "Campylobacter",
#                       species = "jejuni",
                       NCBIFilesDir=".",
                       databaseOnly=TRUE)
## This example demonstrates the limitations of using a lookup table (IOW unsupported characters can end up in the package name) 

## Lets see if this will guess the genus and species here
makeOrgPackageFromNCBI(version="0.1",
                       maintainer="Pengfei Liu <liupfskygre@gmail.com>",
                       author="Pengfei Liu <liupfskygre@gmail.com>",
                       outputDir=".",
                       tax_id="1041930",
#                       genus="Methanocella",
#                       species="conradii",
                       NCBIFilesDir=".")
## And this simply fails to find a match for this species...
## There is no helping this situation since I have not thrown out any of the available scientific names.  Users will simply have to provide a name in this case.



## Barley works...
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "112509",
                       genus = "Hordeum",
                       species = "vulgare",
                       NCBIFilesDir=".")

## Axolotl works...
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "8296",
                       genus = "Ambystoma",
                       species = "mexicanum",
                       NCBIFilesDir=".")






###############################################################################
## Next up: code for telling which things are viable:
## Code should 1) see if NCBI.sqlite 'cache db' exists (& if not -
## make one) and then 2) use the data in there to pre-compute which
## things are viable for making into DBs.
## (Also I need to work out how to know which things I can get data
## from blast2GO for) - pre-download this too?



###############################################################################
## AND refactor code to support the "NOSCHEMA.DB" : done

## Compare DB generated to the one made by makeOrgPackage() -
## NOSCHEMA.DB is better design.

## the older ORGANISM.db schema was unecessarily complicated. - I
## should really upgrade all these things to the new stuff right
## now. - and ALSO: the new schema allows us to make compatible with
## NOSCHEMACHIP.db packages (simpler and more extensible).
## ALSO: the new schema is less complex and better organized since it
## is generated by a hueristic (instead of being bespoke).

## What I need to do here is to support the NOSCHEMA.DB schema and
## still support making the ORGANISM.DB packages via an argument to
## makeOrgPackageFromNCBI.  Old school packages should be born in a
## 'deprecated' state.  (warn people to not use them and the argument
## should also be documented in a way that discourages it's use)

















##########################################################################
## ALSO: While I am in here I need to find a better way to deal with this:
## Basically I need an exception in place for the select code and also
## a more general solution for duplicated stuff
## goids <- keys(mouse4302.db, "GO")
## xx = select(mouse4302.db, head(goids), "GOALL", "GO")
## which takes forever and returns many duplicates
## The above is a select bug that I need to fix in AnnotationDbi...



