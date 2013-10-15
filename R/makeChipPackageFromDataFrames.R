
## function to put together the database.
## This takes a named list of data.frames.
makeChipDbFromDataFrame <- function(probeFrame, tax_id, genus, species, 
                                    dbFileName){
    ## set up DB connection 
    require(RSQLite)
    if(file.exists(dbFileName)){ file.remove(dbFileName) }
    con <- dbConnect(SQLite(), dbFileName)
    AnnotationForge:::.createMetadataTables(con)
    
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


## TODO: Add args for manufacturer information (url, manuf. &
## chipname) - Then pass these along to the seed

## function to make the package:
makeChipPackage <- function(probeFrame, ## data.frame with probe 2 gene mappings
                            orgPkgName, ## name of org package for dependency
                            version,
                            maintainer,
                            author,
                            outputDir = ".",
                            tax_id,
                            genus,
                            species){

    ## probeFrame has to be a data.frame with two columns. (for probes
    ## and gene id)
    
    ## the internal table name should be called probes with fields
    ## (PROBE and GID) to go with new NOSCHEMA style of database.

    ## What if someone wants to use an old style org package with a
    ## new style of chip package?  - That situation could have gotten
    ## complicated BUT chip packages DO specify the org package that
    ## they are supposed to depend on.

    ## SO: when I go to join later on, I will have to check the
    ## schema, and act accordingly by using a helper function. - but
    ## its not that simply since the select method can't reasonably
    ## cross-pollinate

    ## OR: I may need to just require that they use an org package
    ## that is 1) installed and 2) of the NOSCHEMA type.

    ## OR (best idea I think) I may require that the org package be of
    ## the correct type and then make a chip package of the correct
    ## type.

    
    
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
    
    ## generate name from the genus and species
    dbName <- .generateOrgDbName(genus,species)
    ## this becomes the file name for the DB
    dbFileName <- file.path(outputDir,paste0(dbName, ".sqlite"))
    ## Then make the DB
    makeChipDbFromDataFrame(probeFrame, tax_id, genus, species, dbFileName)

    ## make the seed
    seed <- new("AnnDbPkgSeed",
                Package= paste0(dbName,".db"),
                Version=version,
                Author=author,
                Maintainer=maintainer,
                PkgTemplate="NOSCHEMA.DB",  ## TODO: make NOSCHEMACHIP.DB
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
}





