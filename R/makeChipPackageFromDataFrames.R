## For when we have a new schema
.makeProbesTable <- function(probeFrame, con){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS probes (
      PROBEID VARCHAR(80),           -- PROBEID
      GID VARCHAR(10) NULL,           -- GENEID
      is_multiple SMALLINT NOT NULL           -- matches multiple genes?
    );")
  dbGetQuery(con, sql)

  values <- data.frame(probeFrame) 
  sql <- paste("INSERT INTO probes(PROBEID, GID, is_multiple)",
               "VALUES(?,?,?);")
  dbBegin(con)
  dbGetPreparedQuery(con, sql, values)
  dbCommit(con)
  dbGetQuery(con,
                 "CREATE INDEX IF NOT EXISTS Fgenes ON probes (GID)")
  dbGetQuery(con,
                 "CREATE INDEX IF NOT EXISTS Fprobes ON probes (PROBEID)")
  message("probes table filled and indexed")
}

## helper to safely rename things
.tweakMapName <- function(df, initial, replacement){
    if(any(df[["map_name"]]==initial)){
        df[df$map_name==initial,1] <- replacement
        return(df)
    }else{
        return(df)
    }
}

## helper to populate missing map_metadata from org packages for olde templates
.cloneMapMetadata <- function(con, orgPkgName){
    ## 1st we need to extract the map_metadata
    orgCon <- AnnotationDbi:::dbconn(eval(parse(text=orgPkgName)))
    mapValues <- dbGetQuery(orgCon, "SELECT * FROM map_metadata")
    ## then get value for accnum and append it
    orgSchema <- .getOrgSchema(orgPkgName)
    if(orgSchema == "ARABIDOPSIS_DB"){
        accnums <- c("ACCNUM","NA","NA","NA")  ## TODO: these could be passed in
        mapValues <- rbind(mapValues, accNums)
    }else{
        accNums <- mapValues[mapValues$map_name== "ENTREZID",]
        accNums[1] <- "ACCNUM"
        mapValues <- rbind(mapValues, accNums)
        mapValues <- .tweakMapName(mapValues, "ALIAS2EG", "ALIAS2PROBE")
        mapValues <- .tweakMapName(mapValues, "PMID2EG", "PMID2PROBE")
        mapValues <- .tweakMapName(mapValues, "GO2EG", "GO2PROBE")
        mapValues <- .tweakMapName(mapValues, "GO2ALLEGS", "GO2ALLPROBES")
        mapValues <- .tweakMapName(mapValues, "PATH2EG", "PATH2PROBE")
        mapValues <- .tweakMapName(mapValues, "ENSEMBL2EG", "ENSEMBL2PROBE")
    }
    message("Populating map_metadata table:")
    sql<- paste("CREATE TABLE IF NOT EXISTS map_metadata (
      map_name VARCHAR(80) NOT NULL,
      source_name VARCHAR(80) NOT NULL,
      source_url VARCHAR(255) NOT NULL,
      source_date VARCHAR(20) NOT NULL
    );")
    dbGetQuery(con, sql)
    sql <- paste("INSERT INTO map_metadata(map_name, source_name, source_url,
                  source_date)",
                 "VALUES(?,?,?,?);")
    dbBegin(con)
    dbGetPreparedQuery(con, sql, mapValues)
    dbCommit(con)
}


## For maintaining the OLDE style schema
.makeLegacyProbesTable <- function(probeFrame, con, orgPkgName,
                                   accessionsFrame){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS probes (
      probe_id VARCHAR(80),           -- probe ID
      gene_id VARCHAR(10) NULL,           -- Gene ID
      is_multiple SMALLINT NOT NULL           -- matches multiple genes?
    );")
  dbGetQuery(con, sql)

  values <- data.frame(probeFrame) 
  psql <- paste("INSERT INTO probes(probe_id, gene_id, is_multiple)",
               "VALUES(?,?,?);")
  dbBegin(con)
  dbGetPreparedQuery(con, psql, values)
  dbCommit(con)
  dbGetQuery(con,
                 "CREATE INDEX IF NOT EXISTS Fgenes ON probes (gene_id)")
  dbGetQuery(con,
                 "CREATE INDEX IF NOT EXISTS Fprobes ON probes (probe_id)")
  ## now set up correct map_metadata
  .cloneMapMetadata(con, orgPkgName)
  message("probes table filled and indexed")

  ###########################################################
  ## code for legacy accessions table
  dbGetQuery(con, "CREATE TABLE accessions (probe_id VARCHAR(80),accession VARCHAR(20))")
  ## accessionsFrame will either be null or we will need to insert it
  if(!is.null(accessionsFrame)){
      message("Populating accessions table:")
      values <- data.frame(accessionsFrame) 
      asql <- paste("INSERT INTO accessions (probe_id, accession)",
                   "VALUES(?,?);")
      dbBegin(con)
      dbGetPreparedQuery(con, asql, values)
      dbCommit(con)
      ## If there are probes that were mentioned in
      ## accessionsFrame that were NOT mentioned in
      ## probeFrame, that we need to add them like this:
      accProbes <- unique(accessionsFrame[,1])
      pfProbes <- unique(probeFrame[,1])
      if(any(!(accProbes %in% pfProbes))){
          newProbes <- accProbes[!(accProbes %in% pfProbes)]
          newVals <- data.frame(probe_id=newProbes,
                                is_multiple=rep(0L, times=length(newProbes)))
          ssql <- paste("INSERT INTO probes(probe_id, is_multiple)",
                        "VALUES(:probe_id,:is_multiple);")
          ## may have to switch NA to NULL
          dbBegin(con)
          dbGetPreparedQuery(con, ssql, newVals)
          dbCommit(con)
      }
  }
  dbGetQuery(con, "CREATE INDEX Fgbprobes ON accessions (probe_id)")
}


## We need a listing of supported NCBI legacy types
supportedNCBItypes <- function(){
    c("ARABIDOPSIS_DB","ANOPHELES_DB","BOVINE_DB","CANINE_DB","CHICKEN_DB",
      "CHIMP_DB","COELICOLOR_DB","ECOLI_DB","FLY_DB","HUMAN_DB","MALARIA_DB",
      "MOUSE_DB","PIG_DB","RAT_DB","RHESUS_DB","WORM_DB","XENOPUS_DB",
      "ZEBRAFISH_DB")
}


## helper to determine the org schema from the package name string
.getOrgSchema <- function(orgPkgName){
    orgPkg <- eval(parse(text=orgPkgName))
    metadata(orgPkg)[metadata(orgPkg)$name=="DBSCHEMA",][,2]
}
## helper to get the right chip schema for an org schema string
.getChipSchema <- function(orgSchema){
    if(orgSchema == 'NOSCHEMA_DB'){
        chipSchema <- "NOCHIPSCHEMA_DB"
    }else{
        chipSchema <- sub("_DB$","CHIP_DB",orgSchema)
    }
    chipSchema
}
## central ID 
.getChipCentralID <- function(orgPkgName){
    orgPkg <- eval(parse(text=orgPkgName))
    centralID <- metadata(orgPkg)[metadata(orgPkg)$name=="CENTRALID",][,2]
    switch(centralID,
           "EG"="ENTREZID",
           "ORF"="ORF",
           "TAIR"="TAIR",
           "GID"="GID",
           "GID")
}


## Do some work to see which things in probeFrame are multiples
.testForMultiples <- function(probe,probes){
    if(length(probes[probes %in% probe]) > 1){
        1
    }else{
        0
    }
}


## function to put together the database.
## This takes a named list of data.frames.
makeChipDbFromDataFrame <- function(probeFrame, orgPkgName, tax_id,
                                    genus, species, dbFileName,
                                    optionalAccessionsFrame){

    ## 1st connect to the org package
    require(orgPkgName, character.only = TRUE)
    ## Then make final column for the probeFrame            
    multiple <- unlist(lapply(as.character(probeFrame$probes),
                              .testForMultiples,
                              probes=as.character(probeFrame$probes)))
    probeFrame <- cbind(probeFrame, multiple)

    ## set up DB connection 
    require(RSQLite)
    if(file.exists(dbFileName)){ file.remove(dbFileName) }
    con <- dbConnect(SQLite(), dbFileName)
    AnnotationForge:::.createMetadataTables(con)
    
    ## then find out what kind of DB we are building and make matching thing.
    orgSchema <- .getOrgSchema(orgPkgName)
    chipSchema <- .getChipSchema(orgSchema)
    if(orgSchema == 'NOSCHEMA_DB'){
        ## Make probes table
        AnnotationForge:::.makeProbesTable(probeFrame, con)
        ## Add metadata with new schema
        AnnotationForge:::.addEssentialMetadata(con, tax_id, genus, species,
                                                schema=chipSchema,
                                                type="ChipDb",centralID="GID")
    }else if(orgSchema %in% supportedNCBItypes()){
        ## supports the classic "ENTREZID" based legacy org packages
        AnnotationForge:::.makeLegacyProbesTable(probeFrame, con, orgPkgName,
                                                 optionalAccessionsFrame)
        centralID <- .getChipCentralID(orgPkgName)
        AnnotationForge:::.addEssentialMetadata(con, tax_id, genus, species,
                                                schema=chipSchema,
                                                type="ChipDb",
                                                centralID=centralID)
    }else{## note that "legacy" yeast is not supported.
        stop("The org package you have specified has an unsupported schema.")
    }
    ## Add org package info to metadata
    shortPkgName <- sub(".db","",orgPkgName)
    AnnotationForge:::.addMeta(con, "ORGPKGDEP", shortPkgName)
    
}




## function to make the package:
makeChipPackage <- function(prefix,
                            probeFrame, ## data.frame with probe 2 gene mappings
                            orgPkgName, ## name of org package for dependency
                            version,
                            maintainer,
                            author,
                            outputDir = ".",
                            tax_id,
                            genus,
                            species,
                            optionalAccessionsFrame=NULL){

    ## probeFrame has to be a data.frame with two columns. (for probes
    ## and gene id)  And probeFrame HAS to have unique rows
    if(dim(probeFrame)[1] != dim(unique(probeFrame))[1])
        stop("All rows in 'probeFrame' must be unique")
    if(dim(probeFrame)[1] == 0 ||
           dim(probeFrame)[2] != 2){
        stop("'probeFrame' should have two columns with row data for which probes match up with which gene Ids")
        }
    ## check other arguments
    if(!.isSingleString(prefix))
        stop("'prefix' must be a single string")
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
    ## optionalAccessionsFrame
    if(!is.null(optionalAccessionsFrame)){
        if(dim(optionalAccessionsFrame)[1] == 0 ||
           dim(optionalAccessionsFrame)[2] != 2){
            stop("When provided, 'optionalAccessionsFrame' should have two columns with data for which probes match up with which accessions")
        }
        if(dim(optionalAccessionsFrame)[1] != dim(unique(optionalAccessionsFrame))[1])
            stop("All rows in 'optionalAccessionsFrame' must be unique")  
    }
    
    ## The file name for the DB
    dbFileName <- file.path(outputDir,paste0(prefix, ".sqlite"))
    ## Then make the DB
    makeChipDbFromDataFrame(probeFrame, orgPkgName, tax_id,
                            genus, species, dbFileName, optionalAccessionsFrame)

    ## choose the appropriate pkgTemplate (schema)
    orgSchema <- .getOrgSchema(orgPkgName)
    chipSchema <- .getChipSchema(orgSchema)
    ## set up appropriate template 
    if(chipSchema %in% c("NOCHIPSCHEMA_DB", "ARABIDOPSISCHIP_DB")){
        pkgTemplate <- sub("_",".",chipSchema)
    }else{
        pkgTemplate <- "NCBICHIP.DB"
    }
        
    ## make the seed
    seed <- new("AnnDbPkgSeed",
                Package= paste0(prefix,".db"),
                Version=version,
                Author=author,
                Maintainer=maintainer,
                PkgTemplate=pkgTemplate,  
                AnnObjPrefix=prefix,
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
    file.path(outputDir,paste0(prefix,".db"))
}




