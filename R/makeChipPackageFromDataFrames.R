## For when we have a new schema
.makeProbesTable <- function(probeFrame, con){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS probes (
      PROBEID VARCHAR(80),           -- PROBEID
      GENEID VARCHAR(10) NULL,           -- GENEID
      is_multiple SMALLINT NOT NULL           -- matches multiple genes?
    );")
  sqliteQuickSQL(con, sql)

  values <- data.frame(probeFrame) 
  sql <- paste("INSERT INTO probes(PROBEID, GENEID, is_multiple)",
               "VALUES(?,?,?);")
  dbBeginTransaction(con)
  dbGetPreparedQuery(con, sql, values)
  dbCommit(con)
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS Fgenes ON probes (GENEID)")
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS Fprobes ON probes (PROBEID)")
  message("probes table filled and indexed")
}


## helper to populate missing map_metadata from org packages for olde templates
.cloneMapMetadata <- function(con, orgPkgName){
    ## 1st we need to extract the map_metadata
    orgCon <- AnnotationDbi:::dbConn(eval(parse(text=orgPkgName)))
    mapValues <- sqliteQuickSQL(orgCon, "SELECT * FROM map_metadata")
    ## then get value for accnum and append it
    orgSchema <- .getOrgSchema(orgPkgName)
    if(orgSchema == "ARABIDOPSIS_DB"){
        accnums <- c("ACCNUM","NA","NA","NA")  ## TODO: these could be passed in
        mapValues <- rbind(mapValues, accNums)
    }else{
        accNums <- mapValues[mapValues$map_name== "ENTREZID",]
        accNums[1] <- "ACCNUM"
        ## alias2p <- mapValues[mapValues$map_name== "ENTREZID",]
        ## accNums[1] <- "ALIAS2PROBE"
        mapValues <- rbind(mapValues, accNums)
        mapValues[mapValues$map_name== "ALIAS2EG",1] <- "ALIAS2PROBE"
        mapValues[mapValues$map_name== "PMID2EG",1] <- "PMID2PROBE"
        mapValues[mapValues$map_name== "GO2EG",1] <- "GO2PROBE"
        mapValues[mapValues$map_name== "GO2ALLEGS",1] <- "GO2ALLPROBES"
        mapValues[mapValues$map_name== "PATH2EG",1] <- "PATH2PROBE"
        mapValues[mapValues$map_name== "ENSEMBL2EG",1] <- "ENSEMBL2PROBE"
    }
    message("Populating map_metadata table:")
    sql<- paste("CREATE TABLE IF NOT EXISTS map_metadata (
      map_name VARCHAR(80) NOT NULL,
      source_name VARCHAR(80) NOT NULL,
      source_url VARCHAR(255) NOT NULL,
      source_date VARCHAR(20) NOT NULL
    );")
    sqliteQuickSQL(con, sql)
    sql <- paste("INSERT INTO map_metadata(map_name, source_name, source_url,
                  source_date)",
                 "VALUES(?,?,?,?);")
    dbBeginTransaction(con)
    dbGetPreparedQuery(con, sql, mapValues)
    dbCommit(con)
}


## For maintaining the OLDE style schema
.makeLegacyProbesTable <- function(probeFrame, con, orgPkgName){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS probes (
      probe_id VARCHAR(80),           -- probe ID
      gene_id VARCHAR(10) NULL,           -- Gene ID
      is_multiple SMALLINT NOT NULL           -- matches multiple genes?
    );")
  sqliteQuickSQL(con, sql)

  values <- data.frame(probeFrame) 
  sql <- paste("INSERT INTO probes(probe_id, gene_id, is_multiple)",
               "VALUES(?,?,?);")
  dbBeginTransaction(con)
  dbGetPreparedQuery(con, sql, values)
  dbCommit(con)
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS Fgenes ON probes (gene_id)")
  sqliteQuickSQL(con,
                 "CREATE INDEX IF NOT EXISTS Fprobes ON probes (probe_id)")
  ## now set up correct map_metadata
  .cloneMapMetadata(con, orgPkgName)
  
  message("probes table filled and indexed")
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




## function to put together the database.
## This takes a named list of data.frames.
makeChipDbFromDataFrame <- function(probeFrame, orgPkgName, tax_id,
                                    genus, species, dbFileName){

    ## 1st connect to the org package
    require(orgPkgName, character.only = TRUE)
    
            
    ## Do some work to see which things in probeFrame are multiples
    testForMultiples <- function(probe,probes){
        if(length(probe %in% probes) > 1){
            1
        }else{
            0
        }
    }
    multiple <- unlist(lapply(as.character(probeFrame$probes),
                              testForMultiples,
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
        AnnotationForge:::.makeLegacyProbesTable(probeFrame, con, orgPkgName)
        centralID <- .getChipCentralID(orgPkgName)
        AnnotationForge:::.addEssentialMetadata(con, tax_id, genus, species,
                                                schema=chipSchema,
                                                type="ChipDb",
                                                centralID=centralID)
    }else{## note that "legacy" yeast is not supported.
        stop("The org package you have specified has an unsupported schema.")
    }
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
                            species){

    ## probeFrame has to be a data.frame with two columns. (for probes
    ## and gene id)  And probeFrame HAS to have unique rows
    if(dim(probeFrame)[1] != dim(unique(probeFrame))[1])
        stop("All rows in 'probeFrame' must be unique")
    
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
    
    ## The file name for the DB
    dbFileName <- file.path(outputDir,paste0(prefix, ".sqlite"))
    ## Then make the DB
    makeChipDbFromDataFrame(probeFrame, orgPkgName, tax_id,
                            genus, species, dbFileName)

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




## replacing map_metadata worked!  But stuff is STILL missing.
## TODO: still missing: @ALIAS2PROBESOURCE@
