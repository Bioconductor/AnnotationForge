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

