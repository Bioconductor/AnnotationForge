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
## select().  Perhaps NO_SCHEMA? (to indicate no predefined schema)

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

## OR, it could be deduced if the user is required to name their
## passed in data frames the same as the thing that they want back out
## at the end.

## SO IOW a table like this (which also has to itself be named in the list):

### GENE_ID  SYMBOL
### 1        msx2
### 2        hoxa13

## Would mean that cols would return SYMBOL.  ALSO, we have to
## standardize on the name of "GENE_ID" (or on 1st column) to be the
## way the df indicates the primary key for the DB.

## An internally, that table would look like this:

### _id      GENE_ID
### 1        1
### 2        2

### _id      SYMBOL
### 1        msx2
### 2        hoxa13

## And the 1st table is special.  It just "grows" as more unique
## GENE_IDs show up

## This would all work with the one exception being that if they give
## me GO, they have to also give me EVIDENCE and ONTOLOGY (something).
## So my code would have to check for that (IOW, if GO is detected,
## there should also be a field called ONTOLOGY and one called
## EVIDENCE, and they should share a table) - but I don't think I
## should check for that?  Or now I have to since I return GO terms
## expecting to find those fields as well.

## TODO:
## OR: make the code that gets those GO terms more permissive (IOW
## just "try" to get the extra fields - but don't require it).





## This makes the special genes table of an EG DB.
.makeGenesTable <- function(genes, con){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS genes (
      _id INTEGER PRIMARY KEY,
      GID VARCHAR(10) NOT NULL UNIQUE           -- Gene ID
    );")
  sqliteQuickSQL(con, sql)

  geneid <- data.frame(genes) ## TODO: data.frame() necessary???
  sql<- paste("INSERT INTO genes(GID) VALUES(?);")
  dbBeginTransaction(con)
  dbGetPreparedQuery(con, sql, geneid)
  dbCommit(con)
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS genes__id_ind ON genes (_id)")    
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS genes_GID_ind ON genes (GID)")
  message("genes table filled")
}



.makeTable <- function(data, table, con, fieldNameLens=25){
    ## indFields tracks the things to index (not GID but _id, plus the rest)
    indFields <- c(names(data)[!(names(data) %in% "GID")],"_id")
    message(paste("Populating",table,"table:"))
    tableFieldLines <- paste(paste(names(data)[-1]," VARCHAR(",
                                 fieldNameLens,") NOT NULL,    -- data"),
                           collapse="\n       ")
  ## For temp table, lets do it like this:
  if(dim(data)[1] == 0){
    ## if we don't have anything to put into the table, then we don't even
    ## want to make a table.
    warning(paste("no values found for table ",table,
                  " in this data chunk.", sep=""))
    ## Create our real table.
    .makeEmptySimpleTable(con, table, tableFieldLines)
    return()
  }else{
    dbWriteTable(con, "temp", data, row.names=FALSE)
    ## Create our real table.
    .makeEmptySimpleTable(con, table, tableFieldLines)
    selFieldLines <- paste(paste("t.",names(data)[-1],sep=""),collapse=",")
    sql<- paste0("
    INSERT INTO ",table,"
     SELECT g._id as _id, ",selFieldLines,"
     FROM genes AS g, temp AS t
     WHERE g.GID=t.GID
     ORDER BY g._id;")
    sqliteQuickSQL(con, sql)

    ## Add index to all fields in indFields (default is all)
    for(i in seq_len(length(indFields))){
    sqliteQuickSQL(con,
        paste0("CREATE INDEX IF NOT EXISTS ",
              table,"_",indFields[i],"_ind ON ",table,
              " (",indFields[i],");"))    
    }
    
    ## drop the temp table
    sqliteQuickSQL(con, "DROP TABLE temp;")
  }
  message(paste(table,"table filled"))
}



## helper to add most basic of metadata
.addEssentialMetadata <- function(con, tax_id, genus, species){
  name <- c("DBSCHEMAVERSION","DBSCHEMA","ORGANISM","SPECIES","CENTRALID",
            "TAXID",
            "Db type","Supporting package")
  value<- c("2.1","NOSCHEMA_DB",paste(genus,species),paste(genus,species),
            "GID",tax_id,
            "OrgDb","AnnotationDbi")
  AnnotationForge:::.addMeta(con, name, value)
}

## helper to prepare/filter data for two GO tables.
.makeNewGOTables <- function(con, goTable, goData){
    ## So 1st drop the old go table
    sqliteQuickSQL(con, paste0("DROP TABLE ",goTable,";"))
    ## And then make the new ones.
    require(GO.db)
    goOnts <- select(GO.db, as.character(goData$GO), "ONTOLOGY")
    if(dim(goOnts)[1] == dim(goData)[1]){
        goData <- cbind(goData, goOnts)
    }else{stop("Ontology should be 1:1 with GOIDs")}
    ## and re-shuffle
    goData <- goData[,c("GID","GO","EVIDENCE","ONTOLOGY")]
    ## now filter that for terms that are "too new"
    message("Dropping GO IDs that are too new for the current GO.db")
    goData <- goData[goData[["GO"]] %in% Lkeys(GOTERM),]
    ## Then make the 1st table
    .makeTable(goData, names(goData), con=con)
     
    
    
}


## function to put together the database.
## This takes a named list of data.frames.
makeOrgDbFromDataFrames <- function(data, tax_id, genus, species, 
                                    dbFileName, goTable){    
    ## set up DB connection 
    require(RSQLite)
    if(file.exists(dbFileName)){ file.remove(dbFileName) }
    con <- dbConnect(SQLite(), dbFileName)
    AnnotationForge:::.createMetadataTables(con)
    ## TODO: why can't I drop these (investigate)
#    sqliteQuickSQL(con, "DROP TABLE map_counts;")
#    sqliteQuickSQL(con, "DROP TABLE map_metadata;")
    
    ## gather all GIDs together and make the genes table
    genes <- unique(unlist(unname(lapply(data, "[", 'GID'))))
    AnnotationForge:::.makeGenesTable(genes, con)
    
    ## Then do each data.frame in turn
    mapply(FUN=.makeTable, data, names(data), MoreArgs=list(con=con))
    
    
    
    ## Add metadata but keep it very basic
    AnnotationForge:::.addEssentialMetadata(con, tax_id, genus, species)
        
    ## when we have a goTable, we make special GO tables
    if(goTable %in% names(data)){
        ## An extra check for go table (when specified)
        goData <- data[[goTable]]
        if(!all(names(goData) == c("GID", "GO", "EVIDENCE")))
      stop("'goTable' must have three columns called 'GID','GO' and 'EVIDENCE'")
        .makeNewGOTables(con, goTable, goData)
    }

}


## TODO: change the function so it takes ... instead of list (so that
## the arguments are named).  This can be switched from the list later
## by just trapping what is in ... and passing it on as a named list.
## like this:  data <- list(...)

## TODO: this should return the path to the created dir so that the
## user can just call install.packages(on that filepath).

.isSingleString <- function(x){
    is.atomic(x) && length(x) == 1L && is.character(x)
}
.isSingleStringOrNA <- function(x)
{
    is.atomic(x) && length(x) == 1L && (is.character(x) || is.na(x))
}


## function to make the package:
makeOrgPackage <- function(data,
                           version,
                           maintainer,
                           author,
                           outputDir = ".",
                           tax_id,
                           genus,
                           species,
                           goTable=NA){

    ## Data has to meet some strict criteria.
    ## check that it's a list of data.frames.
    dataClasses <- unique(sapply(data, class))
    if(!is.list(data) || dataClasses!="data.frame")
        stop("'data' must be a named list of data.frame objects")
    ## There must be unique names for each data.frame. (
    if(length(unique(names(data))) != length(data))
        stop("All elements of 'data' must be a named")
    ## None of the list names is allowed to be "genes", "metadata" 
    blackListedNames <- c("genes","metadata")
    if(any(names(data) %in% blackListedNames))
       stop("'genes' and 'metadata' are reserved.  Please choose different names for elements of 'data'.")
    ## The data.frames should NOT be allowed to have missing/redundant rows...
    lengthsUni <- sapply(data, function(x){dim(unique(x))[1]})
    lengthsRaw <- sapply(data, function(x){dim(x)[1]})
    if(any(lengthsUni != lengthsRaw))
        stop("'data' should not contain redundant rows")
    ## There must be colnames for each data.frame.(they must each be
    ## unique)
    .noGID <- function(x){x[!(x %in% "GID")]}
    colnamesUni <- unique(.noGID(unlist(sapply(data, colnames))))
    colnamesAll <- .noGID(unlist(sapply(data, colnames)))
    names(colnamesAll) <- NULL
    if(any(colnamesUni != colnamesAll))
        stop("data.frames in 'data' should have unique names for all fields that are not the primary gene id 'GID'")
    ## The 1st column of each data.frame must be a gene ID  (GID)
    colnameGIDs <- sapply(data, function(x){colnames(x)[1]})
    if(any(colnameGIDs != "GID"))
        stop("The 1st column must always be the gene ID 'GID'")
    ## The GID columns all have to be the same type.
    GIDCols <- unique(sapply(data, function(x){class(x[["GID"]])}))
    if(length(GIDCols) >1)
        stop("The type of data in the 'GID' columns must be the same for all data.frames")
    ## check other arguments
    if(!.isSingleString(version))
        stop("'version' must be a single string")
    if(!.isSingleString(maintainer))
        stop("'maintainer' must be a single string")
    if(!.isSingleString(author))
        stop("'author' must be a single string")
    if(outputDir!="." && file.access(outputDir)[[1]]!=0){
        stop("Selected outputDir '", outputDir,"' does not exist.")}
    if(!.isSingleString(tax_id))
        stop("'tax_id' must be a single string")
    if(!.isSingleString(genus))
        stop("'genus' must be a single string")
    if(!.isSingleString(species))
        stop("'species' must be a single string")
    if(!.isSingleStringOrNA(goTable))
        stop("'goTable' argument needs to be a single string or NA")
    if(!is.na(goTable) && !(goTable %in% names(data)))
        stop("When definined, 'goTable' needs to be a table name from 'data'")
    
    
    ## generate name from the genus and species
    dbName <- .generateOrgDbName(genus,species)
    ## this becomes the file name for the DB
    dbFileName <- file.path(outputDir,paste0(dbName, ".sqlite"))
    ## Then make the DB
    makeOrgDbFromDataFrames(data, tax_id, genus, species, dbFileName, goTable)
        
  
    seed <- new("AnnDbPkgSeed",
                Package= paste0(dbName,".db"),
                Version=version,
                Author=author,
                Maintainer=maintainer,
                PkgTemplate="NOSCHEMA.DB",
                AnnObjPrefix=dbName,
                organism = paste(genus, species),
                species = paste(genus, species),
                biocViews = "annotation",
                manufacturerUrl = "no manufacturer",
                manufacturer = "no manufacturer",
                chipName = "no manufacturer")
    
    makeAnnDbPkg(seed, dbFileName, dest_dir=outputDir)
    
    ## cleanup
    message("Now deleting temporary database file")
    file.remove(dbfile)
    ## return the path to the dir that was just created.
    file.path(outputDir,paste0(dbName,".db"))
}





