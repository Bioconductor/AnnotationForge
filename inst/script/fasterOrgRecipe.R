## The aim here is to sketch out and then implement a strategy for
## making all the org packages quickly by calling the makeOrgPackage()
## on code processed from the FTP file downloaded by the NCBI_ftp.R
## source.

## GOALS:
## 1) I want this to be performant.  It cannot take 30 minutes to
## download the whole FTP file and split that every time (and most of
## the time is in the splitting).  I need to do the initial processing
## step once and then run the recipe on some generated stuff that
## comes out of that.
## 2) I want to pre-filter based on a. whether the tax ID is legit,
## b. whether there are records for that id and c. whether there are
## GO annotations for that organism either at NCBI OR at blast2GO.
## 3) It might makes sense to use the makeOrgPackage code as a base
## for making the classic org packages FIRST (it's basically an good
## opportunity to refactor).

## Things to do:
## 1) refactor makeOrgDbFromNCBI() to use makeOrgPackage() internally.
## 2) make another function makeOrgDbFromNCBIFile() that starts with a
## saved file of NCBI data and a saved block of GO data (that each
## contain only one TaxID worth of code).  Modify makeOrgDbFromNCBI()
## so that it uses this function for the 2nd half of what it needs to
## do (after it does a one off of the slow step).
## 3) make a separate helper to be called 1st that just processes
## one/all of the NCBI tax IDs out of the file above (as a saved file
## that can be fed to makeOrgDbFromNCBIFile()) so that we can make one
## or many of these things depending on what is specified.  This is
## the place to also do pre-filtering.  No sense saving what you
## cannot process later.
## 4) modify the recipe so that I can do a pre-process step.  I should
## probably make sure that I can do this FIRST since it is crucial to
## making this all work. - And presently: I cannot do this.  Or I
## *could*, but not without messing up the metadata...

## SO: I might need a different strategy.  I might want to just find a
## better way to pre-compute the things that I will have available as
## NCBI data AND as GO records so that I know that up front... 




## So lets explore an example and see what is actually the slow part...
## slow step is the part where we discard data.


debug(makeOrgDbFromNCBI)
debug(AnnotationForge:::.makeBaseDBFromDLs)
debug(AnnotationForge:::.saveFiles)
debug(AnnotationForge:::.downloadData)

debug(AnnotationForge:::.getFiles)


debug(AnnotationForge:::.writeToNCBIDB)

debug(AnnotationForge:::.indexTaxIds)


library(AnnotationForge)

## use the NCBIFilesDi arg to cache the big file locally
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Some One <so@someplace.org>",
                       maintainer = "Some One <so@someplace.org>",
                       outputDir = ".",
                       tax_id = "192222",
                       genus = "Campylobacter",
                       species = "jejuni",
                       NCBIFilesDir=".")



## What I really want to do is just do more inside of .downloadData()
## so that instead of only saving the data for ONE organism, I instead
## save it for all organisms, but named like "txid.txt".  That way I
## only have to do the preprocess step the 1st time...
## NO: LEAVE it.  Right now it's better for memory while being about the
## same for speed as reading in once.

## OR I might be able to just work from the one file.  SLOW part is
## read.delim call (and I am already doing that as fast as possible.
## Speedup: 1) only call the slow part ONCE per critter (not twice).
## 2) be sure the source file is saved in FULL and not overwritten.
## NO this won't work unless we change the code for AnnotationHubData: a lot.

## To learn which things are viable as tax IDs, I need to just DL all
## the files and do an intersection of the uniques from their Tax_IDs.


## You know now that I look at it more, I think it might be worthwhile
## to save each main file as a table in a core DB called NCBI (which
## could be made if it doesn't exist already in the dir.  Some care
## will have to be taken to remove this dir if it is more than a week
## old (it could be date stamped).  But it could have data in the form
## of tables from each of the core files we download (with indices for
## the Tax_id).  And then getting out the data from disc each time
## could just be a select operation.

