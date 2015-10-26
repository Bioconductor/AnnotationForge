## This makes the special genes table of an EG DB.
.makeGenesTable <- function(genes, con){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS genes (
      _id INTEGER PRIMARY KEY,
      GID VARCHAR(10) NOT NULL UNIQUE           -- Gene ID
    );")
  dbGetQuery(con, sql)

  geneid <- data.frame(genes) ## TODO: data.frame() necessary???
  sql<- paste("INSERT INTO genes(GID) VALUES(?);")
  dbBegin(con)
  dbGetPreparedQuery(con, sql, geneid)
  dbCommit(con)
  dbGetQuery(con,
                 "CREATE INDEX IF NOT EXISTS genes__id_ind ON genes (_id)")    
  dbGetQuery(con,
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
    dbGetQuery(con, sql)

    ## Add index to all fields in indFields (default is all)
    for(i in seq_len(length(indFields))){
    dbGetQuery(con,
        paste0("CREATE INDEX IF NOT EXISTS ",
              table,"_",indFields[i],"_ind ON ",table,
              " (",indFields[i],");"))    
    }
    
    ## drop the temp table
    dbGetQuery(con, "DROP TABLE temp;")
  }
  message(paste(table,"table filled"))
}



## helper to add most basic of metadata
.addEssentialMetadata <- function(con, tax_id, genus, species,
                                  schema="NOSCHEMA_DB",
                                  type="OrgDb",
                                  centralID="GID"){
  name <- c("DBSCHEMAVERSION","DBSCHEMA","ORGANISM","SPECIES","CENTRALID",
            "Taxonomy ID",
            "Db type","Supporting package")
  value<- c("2.1",schema,paste(genus,species),paste(genus,species),
            centralID,tax_id,
            type,"AnnotationDbi")
  AnnotationForge:::.addMeta(con, name, value)
}

.addOntologyData <- function(data){
    ## And then make the new ones.
    require(GO.db)
    goOnts <- select(GO.db, as.character(data$GO), "ONTOLOGY")
    if(dim(goOnts)[1] == dim(data)[1]){
        data <- cbind(data, goOnts)
    }else{stop("Ontology should be 1:1 with GOIDs")}
    ## and re-shuffle
    data[,c("GID","GO","EVIDENCE","ONTOLOGY")]
}

## helper to prepare/filter data for two GO tables.
.makeNewGOTables <- function(con, goTable, goData){
    ## So 1st drop the old go table
    dbGetQuery(con, paste0("DROP TABLE ",goTable,";"))
    ## add Ontologies data
    goData <- .addOntologyData(goData)
    ## now filter that for terms that are "too new"
    message("Dropping GO IDs that are too new for the current GO.db")
    goData <- goData[goData[["GO"]] %in% Lkeys(GOTERM),]
    ## Then make the 1st table
    .makeTable(goData, "go", con=con)

    ## Then prepare data for the 2nd (go_all) table
    gbp <- goData[goData$ONTOLOGY=="BP",c("GID","GO","EVIDENCE")]
    gcc <- goData[goData$ONTOLOGY=="CC",c("GID","GO","EVIDENCE")]
    gmf <- goData[goData$ONTOLOGY=="MF",c("GID","GO","EVIDENCE")]
    names(gbp) <- names(gcc) <- names(gmf) <- c("gene_id","go_id","evidence")
    ## then recycle older expand function
    bpAll <- .expandGOFrame(gbp, GOBPANCESTOR)
    mfAll <- .expandGOFrame(gmf, GOMFANCESTOR)
    ccAll <- .expandGOFrame(gcc, GOCCANCESTOR)
    ## then combine
    goAllData <- rbind(bpAll,ccAll,mfAll)
    names(goAllData) <- c("GID","GO","EVIDENCE")
    goAllData <- .addOntologyData(goAllData)
    names(goAllData) <- c("GID","GOALL","EVIDENCEALL","ONTOLOGYALL")
    ## Then make the 2nd table
    .makeTable(goAllData, "go_all", con=con)
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
#    dbGetQuery(con, "DROP TABLE map_counts;")
#    dbGetQuery(con, "DROP TABLE map_metadata;")
    
    ## gather all GIDs together and make the genes table
    genes <- unique(unlist(unname(lapply(data, "[", 'GID'))))
    AnnotationForge:::.makeGenesTable(genes, con)
    
    ## Then do each data.frame in turn
    mapply(FUN=.makeTable, data, names(data), MoreArgs=list(con=con))
        
    ## Add metadata but keep it very basic
    AnnotationForge:::.addEssentialMetadata(con, tax_id, genus, species)
        
    ## when we have a goTable, we make special GO tables
    if(goTable %in% names(data)){
        ## Extra checks for go table (when specified)
        goData <- data[[goTable]]
        if(!all(names(goData) == c("GID", "GO", "EVIDENCE")))
      stop("'goTable' must have three columns called 'GID','GO' and 'EVIDENCE'")
        if(any(!grepl("^GO:", as.character(goData$GO))))
            stop("'goTable' GO Ids must be formatted like 'GO:XXXXXXX'")
        .makeNewGOTables(con, goTable, goData)
    }

}


## TODO: change the function so it takes ... instead of list (so that
## the arguments are named).  This can be switched from the list later
## by just trapping what is in ... and passing it on as a named list.
## like this:  data <- list(...)

## TODO: this should return the path to the created dir so that the
## user can just call install.packages(on that filepath).

.makeOrgPackage <- function(data,
                            version,
                            maintainer,
                            author,
                            outputDir = getwd(),
                            tax_id,
                            genus=NULL,
                            species=NULL,
                            goTable=NA,
                            databaseOnly=FALSE){

    ## unique names for all 'data'
    if (length(unique(names(data))) != length(data))
        stop("All elements of '...' must have unique names")
    blackListedNames <- c("genes","metadata")
    if (any(names(data) %in% blackListedNames))
       stop("'genes' and 'metadata' are reserved. Please choose different ",
            "names for elements of '...'.")
    ## coerce to data.frame
    data <- lapply(data, as.data.frame)

    ## drop rownames, no duplicated rows
    data <- lapply(data, function(x){
                rownames(x) <- NULL
                if (any(duplicated(x)))
                    stop("data.frames in '...' cannot contain duplicated rows")
                x
            })

    ## unique colnames for each data.frame
    .noGID <- function(x){x[!(x %in% "GID")]}
    colnamesUni <- unique(.noGID(unlist(sapply(data, colnames))))
    colnamesAll <- .noGID(unlist(sapply(data, colnames)))
    names(colnamesAll) <- NULL
    if(any(colnamesUni != colnamesAll))
        stop("data.frames should have completely unique names for all fields that are not the primary gene id 'GID'")
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
    if(!.isSingleStringOrNull(genus))
        stop("'genus' must be a single string or NULL")
    if(!.isSingleStringOrNull(species))
        stop("'species' must be a single string or NULL")
    if(!.isSingleStringOrNA(goTable))
        stop("'goTable' argument needs to be a single string or NULL") ## only an NA internally - a NULL is what would have come in from outside...
    if(!is.na(goTable) && !(goTable %in% names(data)))
        stop("When definined, 'goTable' needs to be a table name from the named data.frames passed in to '...'")
    
    ## if genus or species are null, then we should get them now.
    if(is.null(genus)){genus <- GenomeInfoDb:::.lookupSpeciesFromTaxId(tax_id)[['genus']] }
    if(is.null(species)){
        species <- GenomeInfoDb:::.lookupSpeciesFromTaxId(tax_id)[['species']]
        species <- gsub(' ','.', species)
    }

    
    ## generate name from the genus and species
    dbName <- .generateOrgDbName(genus,species)
    ## this becomes the file name for the DB
    dbFileName <- file.path(outputDir,paste0(dbName, ".sqlite"))
    ## Then make the DB
    makeOrgDbFromDataFrames(data, tax_id, genus, species, dbFileName, goTable)
        
    if(databaseOnly==FALSE){
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
        file.remove(dbFileName)
        ## return the path to the dir that was just created.
        file.path(outputDir,paste0(dbName,".db"))
    }else{
        ## return the path to the database file
        file.path(outputDir,dbFileName)
    }
}



## function to make the package:
makeOrgPackage <- function(...,
                           version,
                           maintainer,
                           author,
                           outputDir = getwd(),
                           tax_id,
                           genus=NULL,
                           species=NULL,
                           goTable=NULL){
    if(is.null(goTable)){goTable <- NA}
    ## get all the arguments into a list
    data <- list(...)
    .makeOrgPackage(data,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    outputDir=outputDir,
                    tax_id=tax_id,
                    genus=genus,
                    species=species,
                    goTable=goTable)
}





