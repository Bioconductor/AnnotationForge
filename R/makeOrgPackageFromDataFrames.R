## Function to make an org package from a named list of data.frame
## objects.  Each data.frame should be processed in turn into a table
## in the database.  Each data.frame should have a key named _id that
## allows them to all be linked together.  The named of the list
## elements will become the table names.

## That part is straightforward.  But what about when the user has a
## new kind of table? They should be able to add it here of course.
## But how will that get exposed for select()?  I think there should
## be code that finds new tables at load time and creates symbols
## dynamically for use with select()

## This will happen by a simple addition to .defTables that will add
## new SYMBOLS in an automagical way (following a simple hueristic).

## What class will these be?  OrgDb.  So that case (in
## .definePossibleTables) needs to "scan" for "new" tables.

## What schema? ORGANISM_DB?  We might want to be even less demanding
## (and make a new one that has almost nothing in it).  The basic idea
## is that we won't have anything in there except support for
## select().

## Once we have all the data in a list of data.frames, we will require
## very little in terms of tables.  Just two tables and one of them
## needs to be a genes table (for a proper ID).

## The gene_info table is a weird exception though since it has
## multiple fields that are accessed independently (unlike most of the
## others). - So do I even try to recreate this?  I think I don't want
## to.  I think I want to just do these things separately.  However
## this table is an exceptional case (b/c one record per _id).  There
## will be tables where the data is destroyed by splitting it up.


## For example what about GO?  GO is stored as three tables in the
## schema and then recombined into a GO view.  So that means for
## select, I only really need a single go table in this schema (and a
## go_all table).  BUT: what about the fact that it has multiple
## fields?  Do I need to replicate that too?  YES I do.  Because you
## can't split these fields apart and then reconstitute them later on.

## So that means I need to support multi-field tables too!  But I can
## do that too.  Because I can just add a new line to the .defTables
## for each field.

## So I will need to code for select() to check the schema and do
## things a wee bit differently if I have one of these custom
## ORGANISM_DB objects.

## Lets look at one and see what they look like right now. Right now:
## .defTables is defined to be a stock set of values UNLESS you are
## one of the "stockSpecies", This will have to be changed so that
## instead it 1) checks the schema and if it's ORGANISM_DB, then use
## the standard set otherwise, call out to another helper and try to
## extract it from the DB.

## so as we write the DB, we should also write records into a DB table
## that records the defTables relationships.  IOW we need to record
## like: "ENTREZID" = c("genes","gene_id") etc. in a special three col
## table.  This can't be fully deduced from the colnames of the
## data.frame and the names of the tables already passed in, so we
## will need a special data.frame for this.

## This would all work with the one exception being that if they give
## me GO, they have to also give me EVIDENCE and ONTOLOGY (something).
## So my code would have to check for that (IOW, if GO is detected,
## there should also be a field called ONTOLOGY and one called
## EVIDENCE, and they should share a table)


## function to make the package:
makeOrgPackage <- function(data,
                           fields,
                           version,
                           maintainer,
                           author,
                           outputDir = ".",
                           tax_id,
                           genus,
                           species,
                           NCBIFilesDir=NULL){

  if(outputDir!="." && file.access(outputDir)[[1]]!=0){
    stop("Selected outputDir '", outputDir,"' does not exist.")}

  ## 'outputDir' is not passed to makeOrgDbFromNCBI(). Hence the db file
  ## is always created in ".". Maybe that could be revisited.
  makeOrgDbFromDataFrames(data)
  
  dbName <- .generateOrgDbName(genus,species)
  dbfile <- paste(dbName, ".sqlite", sep="")

  seed <- new("AnnDbPkgSeed",
              Package= paste(dbName,".db",sep=""),
              Version=version,
              Author=author,
              Maintainer=maintainer,
              PkgTemplate="ORGANISM.DB",
              AnnObjPrefix=dbName,
              organism = paste(genus, species),
              species = paste(genus, species),
              biocViews = "annotation",
              manufacturerUrl = "no manufacturer",
              manufacturer = "no manufacturer",
              chipName = "no manufacturer")
  
  makeAnnDbPkg(seed, dbfile, dest_dir=outputDir)
  
  ## cleanup
  file.remove(dbfile)
}
