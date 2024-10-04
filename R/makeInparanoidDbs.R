## This file will contain the code to make ONE Inparanoid DB from
## source files it downloads, unzips and imports (including metadata).
## The code to wrap that file into a package will still live in
## AnnotationForge for now.

## TODO add a table with taxIDs mapped onto the short names.
## this file is here:  http://inparanoid.sbc.su.se/download/current/sequences/species.mapping.inparanoid8

## helper to get a single dir listing.
.getSubDirs <- function(dir)
{
    ## now get the dir
    loadNamespace("XML")
    ## So I need to make a handler that collects the parts I want:
    getLinks = function() {
        links = character()
        list(a = function(node, ...) {
            links <<- c(links, XML::xmlGetAttr(node, "href"))
            node
        },
        links = function()links)
    }
    ## now get the links (sub dirs)
    h1 = getLinks()
    XML::htmlTreeParse(dir, handlers = h1)
    res <- h1$links() ## base result
    ## some filtering
    res <- res[!(res %in% c("?C=N;O=D", "?C=M;O=A", "?C=S;O=A", "?C=D;O=A",
                            "/download/current/"))]
    res
}


## helper to get one file set
.getAFileSet <- function(dir, dataDir){
    files <- .getSubDirs(dir)
    ## filter empty dirs (are not tarballs)
    files <- files[grepl(".tgz$",files)]
    allPaths <- paste0(dir, files)
    ## now actually download all these files
    localFiles <- file.path(dataDir, files)
    mapply(download.file, url=allPaths, destfile=localFiles)
    ## untar all that stuff...
    lapply(localFiles, untar, exdir=dataDir)
    ## And then return the list of junk you just downloaded
    localFiles
}



## helper for removing junk from data and outputting warnings as needed
.prepareData <- function(file, cols) {
    ## read in the data
    ## AND some of their files have a period randomly stuck onto the end...
    ## OKAY then we will do this nonsense to try and recover...
    if(!file.exists(file) && .Platform$OS.type == "unix"){
        message("Attempting to guess the fileName based on:", file)
        file <- system(paste0("ls ",file,"*"), intern=TRUE)
        message("Our guess is that you wanted:", file)
    }
    message("reading in ", file)
    insVals <- read.delim(file=file, header=FALSE, sep="\t", quote="",
                          colClasses=c("integer", rep("character", 5)),
                          stringsAsFactors=FALSE)
    
    ## check to see if any critical stuff is missing
    countCol <- function(col){
        numNA <- sum(is.na(col))
    }
    critVals <- insVals[, cols]
    ## FIXME: don't need to coerce to matrix to get column count of NAs
    ## should be sapply(critVals, countCol)
    NAColCnts <- apply(as.matrix(critVals), 2, countCol)

    ## Then we need to make log entries for flaws that we find...
    clnVals <- insVals
    for(i in 1:length(NAColCnts)){             
        if(NAColCnts[i]>0){
            warning(paste("CRITICAL TABLE FLAW!  There were ",NAColCnts[i],
             " NAs inside of critical col ",cols[i]," inside the file named ",
             file,'\n',sep=""))
          ## cat(paste("CRITICAL TABLE FLAW!  There were ",NAColCnts[i],
          ##  " NAs inside of critical col ",cols[i]," inside the file named ",
          ##  file,'\n',sep=""),file="BADINPSrcFiles.log", append=TRUE)
          ## then scrub out the bad data rows (on crit cols)  
            clnVals <- clnVals[!is.na(insVals[, cols[i]]), ]
        }
    }
    clnVals
}

## helper to get us a 
.lookupTableName <- function(shortName){
    allNames <- read.delim(system.file('extdata','inp8_Full_species_mapping',
                                        package='AnnotationForge'),
                           sep="\t", header=TRUE, stringsAsFactors=FALSE)
    lidx <-grepl(shortName, allNames$inparanoidSpecies)
    res <- allNames[lidx, 'tableNames']
    #gsub(" ", "_", gsub("-","_", gsub("\\.","",res)))
}

## helper to populate one file to become a table in a database
.popInpTable = function(con, file, species, dataDir="."){
    ## cleanup because there is some peculiarity for some of their files.
    file <- sub("\\..tgz","\\.tgz",file)
    ## extract the tableName we need for the table name etc.
    ## fileSpecies <- sub("^InParanoid.+?-","",file, perl=TRUE)
    fileSpecies <- sub("^.*/InParanoid.+?-","",file, perl=TRUE) ## DEBUG
    fileSpecies <- sub(".tgz$","",fileSpecies)
    ## get tableName by translating the fileSpecies
    tableName <- .lookupTableName(fileSpecies)
    ## tableName <- sub("\\.","_",fileSpecies)
    ## tableName <- sub("-","_",tableName)
    ## get the actual extracted src file name that I need.
    srcFile <- file.path(dataDir, paste0("sqltable.",species,"-",fileSpecies))
    
    ##Make a table
    message(paste0("Creating table: ", tableName))
    sql<- paste0("    CREATE TABLE IF NOT EXISTS ", tableName, " (
        clust_id INTEGER NOT NULL,
        clu2 VARCHAR(10) NOT NULL,
        species VARCHAR(15) NOT NULL,
        score VARCHAR(6) NOT NULL,
        inp_id VARCHAR(30) NOT NULL,  
        seed_status CHAR(4));")
    dbGetQuery(con, sql)
    
    ##Populate it with the contents of the filename
    message(cat(paste0("Populating table: ",tableName)))
    clnVals <- .prepareData(file=srcFile)

    sql<- paste0("    INSERT into ", tableName,
       " (clust_id,clu2,species,score,inp_id,seed_status) VALUES
        (?,?,?,?,?,?)")
    dbBegin(con)
    dbGetQuery(con, sql, unclass(unname(clnVals)))
    dbCommit(con)
}

## Helper to toss out all the generated files...
.cleanupFiles <- function(files, dataDir="."){    
    ## unlink allows wildcards
    fileBase <- sub(file.path(dataDir,""),"",
                    sub("InParanoid","",
                        sub(".tgz","",files)))
    unlink(file.path(dataDir,paste0("*",fileBase,"*")))
    unlink(file.path(dataDir,"*.stdout"))
}

## .makeMapCounts <- function(con, species){
##     tabs <- ## lookup function for "Full_species_mapping"
##     ## Then generate sql inserts
##     sqls <- paste0('INSERT INTO map_metadata ', tabs)    
## }


## TODO: add Organism and species to the metadata!
## helper for filling in the metadtaa table
.makeBasicMetadata <- function(con, species){
    message(paste0("Creating metadata table"))
    sql<- paste0("    CREATE TABLE if not exists metadata (
        name VARCHAR(80) PRIMARY KEY,
        value VARCHAR(255));")
    dbGetQuery(con, sql)
    meta <- read.delim(system.file('extdata','inp8_metadata',
                                   package='AnnotationForge'),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
    species <- sub("_", " ", species)
    meta[meta$name=='ORGANISM','value'] <- species
    meta[meta$name=='SPECIES','value'] <- species
    sql<- paste0("INSERT INTO metadata VALUES (:name, :value)")
    dbBegin(con)
    dbGetQuery(con, sql, unclass(unname(meta)))
    dbCommit(con)
}

## function to make a set of tables.
makeInpDb <- function(dir, dataDir="."){
    ## Start by getting all the data we need
    files <- .getAFileSet(dir, dataDir)
    ## temp hack if you don't want to re-DL the files each time    
    ## files = dir('.', pattern='*.tgz')
    ## set up for a database
    ## the connection is for the KIND of DB
    abbrevSpeciesName <- sub("-.*.tgz$","",files[1])
    ## abbrevSpeciesName <- sub("^InParanoid.","",abbrevSpeciesName)
    abbrevSpeciesName <- sub("^.*/InParanoid.","",abbrevSpeciesName) ## DEBUG!
    DBNAME <- paste0("hom.",
                     .lookupTableName(abbrevSpeciesName),
                     ".inp8",
                     ".sqlite")
    con <- dbConnect(SQLite(),dbname=DBNAME)
    
    ## each DB will need a connection
    for(i in seq_along(files)){
        ## Then start making tables...
        .popInpTable(con, file=files[i], species=abbrevSpeciesName,
                     dataDir=dataDir)    
    }

    ## Don't forget the metadata...
    ## .makeMetadata(con, species=.lookupTableName(abbrevSpeciesName))    
    .makeBasicMetadata(con, species=.lookupTableName(abbrevSpeciesName))
    
    ## And end by unlinking all the data we don't need anymore...
    .cleanupFiles(files, dataDir)
    ## and close the connection
    dbDisconnect(con)
    ## return the name of the DB connection for external saving
    DBNAME
}
