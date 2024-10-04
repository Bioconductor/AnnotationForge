#########################################################################
#########################################################################
####         Helper functions for generating a generic DB            ####
#########################################################################
#########################################################################

## This makes the central table of an EG DB.
## FIXME: confirm no ';' separated multiple GID values
.makeCentralTable <-
    function(entrez, con)
{
    message("Populating genes table:")
    sql <- "
        CREATE TABLE IF NOT EXISTS genes (
            _id INTEGER PRIMARY KEY,
            gene_id VARCHAR(10) NOT NULL UNIQUE           -- Entrez Gene ID
        );"
    dbExecute(con, sql)

    gene_id <- data.frame(entrez) ## TODO: data.frame() necessary???
    sql<- "INSERT INTO genes(gene_id) VALUES(?);"
    dbBegin(con)
    dbExecute(con, sql, unclass(unname(gene_id)))
    dbCommit(con)
    message("genes table filled")
}

## The following takes a data.frame and produces a simple table from that.  It
## expects that the 1st column of that data.frame will always be entrez gene
## IDs.  All fields are assumed to be varchars of size equal to the values in
## fieldNameLens.  TODO: The cols in data have to be named and of equal
## length.  indFields is a character vector with the names of fields that we
## want indexed.  By default only _id will be indexed.
.makeEmptySimpleTable <-
    function(con, table, tableFieldLines)
{
    sql<- paste("
        CREATE TABLE IF NOT EXISTS", table, " (
            _id INTEGER NOT NULL,                         -- REFERENCES genes
        ", tableFieldLines, "
        FOREIGN KEY (_id)
        REFERENCES genes (_id));")
    dbExecute(con, sql)
}

.makeSimpleTable <-
    function(data, table, con, fieldNameLens=25, indFields="_id")
{
    message("Populating ", table, " table:")
    tableFieldLines <- paste(paste(names(data)[-1]," VARCHAR(",
                                   fieldNameLens,") NOT NULL,    -- data"),
                             collapse="\n       ")
    ## For temp table, lets do it like this:
    if (dim(data)[1] == 0) {
        ## if we don't have anything to put into the table, then we don't even
        ## want to make a table.
        warning("no values found for table ", table, " in this data chunk.")
        ## Create our real table.
        .makeEmptySimpleTable(con, table, tableFieldLines)
        return()
    } else {
        dbWriteTable(con, "temp", data, row.names=FALSE)
        ## Create our real table.
        .makeEmptySimpleTable(con, table, tableFieldLines)
        selFieldLines <- paste(paste("t.",names(data)[-1],sep=""),collapse=",")
        sql <- paste0("
            INSERT INTO ", table,"
            SELECT g._id AS _id, ", selFieldLines, "
            FROM genes AS g, temp AS t
            WHERE g.gene_id=t.gene_id
            ORDER BY g._id;");
        dbExecute(con, sql)

        ## Add index to all fields in indFields (default is all)
        for(i in seq_len(length(indFields))) {
            dbExecute(con, paste0("
                CREATE INDEX IF NOT EXISTS ", table, "_", indFields[i], "_ind
                ON ", table, " (", indFields[i], ");"))
        }

        ## drop the temp table
        dbExecute(con, "DROP TABLE temp;")
    }
    message(paste(table,"table filled"))
}

## used to gather ancestor nodes for GO terms
.expandGOFrame <-
    function(frame, AncestMap)
{
    if (dim(frame)[1] ==0) {
        res <- data.frame(gene_id=0, go_id=0, evidence=0)
        return(res[FALSE,])
    }
    ## I want to apply through the original frame and call for the ancestor
    ancList <- mget(as.character(frame$go_id), AncestMap, ifnotfound=NA)
    names(ancList) <- frame$gene_id
    eviCodes <- mget(as.character(frame$go_id), AncestMap, ifnotfound=NA)
    names(eviCodes) <- frame$evidence
    expAncList <- unlist2(ancList)
    expEviCodes <- unlist2(eviCodes)
    extraRows <- data.frame(gene_id=names(expAncList), go_id=expAncList,
                            evidence=names(expEviCodes))
    ##remove rows where go_id="all"
    extraRows <- extraRows[extraRows$go_id != "all",]
    unique(rbind(frame,extraRows))
}

## used to make sure that GO IDs that are NOT in the current GO package do not
## end up in our DB
.filterGOFrame <-
    function(frame)
{
    message("Dropping GO IDs that are too new for the current GO.db")
    frame[frame[["go_id"]] %in% Lkeys(GO.db::GOTERM),]
}

## TODO: modify this so that it no longer unwinds lists...
## used to make the 6 custom GO tables
.makeGOTablesFromNCBI <-
    function(con)
{
    bp <- .filterGOFrame(dbGetQuery(con, "
        SELECT DISTINCT gene_id, go_id, evidence
        FROM gene2go
        WHERE category = 'Process'"))

    mf <- .filterGOFrame(dbGetQuery(con, "
        SELECT DISTINCT gene_id, go_id, evidence
        FROM gene2go
        WHERE category = 'Function'"))

    cc <- .filterGOFrame(dbGetQuery(con, "
        SELECT DISTINCT gene_id, go_id, evidence
        FROM gene2go
        WHERE category = 'Component'"))

    headerNames = c("gene_id","go_id","evidence")
    names(bp) <- headerNames
    names(mf) <- headerNames
    names(cc) <- headerNames

    .makeSimpleTable(bp, table = "go_bp", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))

    .makeSimpleTable(mf, table = "go_mf", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))

    .makeSimpleTable(cc, table = "go_cc", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))

    ## Now expand the three data.frames to incl all ancestor terms
    bp_all <- .expandGOFrame(bp, GO.db::GOBPANCESTOR)
    mf_all <- .expandGOFrame(mf, GO.db::GOMFANCESTOR)
    cc_all <- .expandGOFrame(cc, GO.db::GOCCANCESTOR)

    .makeSimpleTable(bp_all, table = "go_bp_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))

    .makeSimpleTable(mf_all, table = "go_mf_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))

    .makeSimpleTable(cc_all, table = "go_cc_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))
}


#########################################################################
#########################################################################
.downloadAndSaveToTemp <-
    function(url, tmp)
{

    if(httr::http_error(url))
        stop("URL '", url, "' does not exist")    

##    binRes <- RCurl::getBinaryURL(url)
##    writeBin(binRes, con=tmp)
    ## RCurl on Windows has TLS problems
    ## Use download.file instead
    ## f = RCurl::CFILE(tmp, mode="wb")
    ##     RCurl::curlPerform(url = url, writedata = f@ref)
    ##     RCurl::close(f)
    options(timeout = max(1000, getOption("timeout")))
    download.file(url, tmp, quiet = TRUE, mode = "wb")
 }

.tryDL <-
    function(url, tmp)
{
    times = 4 ## try this many times to DL
    for(i in seq_len(times)) {
        ## tryResult <- try( download.file(url, tmp, quiet=TRUE) , silent=TRUE)
        tryResult <- try( .downloadAndSaveToTemp(url, tmp) , silent=TRUE)
        if (is(tryResult,"try-error") && i < times) {
            Sys.sleep(20)
        } else if (is(tryResult,"try-error") && i >= times) {
            msg = c("url access failed after ", i, " attempts; url: ", url)
            stop(paste(strwrap(msg, exdent=2), collapse="\n"))
        } else {
            return()
        }
    }
}

.writeToNCBIDB <-
    function(NCBIcon, tableName, filepath, file)
{
    colClasses <- rep("character", times= length(unlist(file)))
    nrows <- 1000000
    overwrite <- TRUE
    append <- FALSE
    ## first chunk
    filecon <- file(description=filepath, open="r")
    vals <- read.delim(filecon, header=FALSE, sep="\t", skip=1, nrows=nrows,
                       stringsAsFactors=FALSE, quote="", colClasses=colClasses)
    repeat {
        if (nrow(vals) == 0)
            break
        ## process
        if ('tax_id' %in% colnames(vals))
            vals[['tax_id']] <- as.integer(vals[['tax_id']])
        colnames(vals) <- unlist(file)
        dbWriteTable(NCBIcon, tableName, vals,
                     overwrite=overwrite, append=append)
        ## break if last chunk was final chunk
        if (nrow(vals) != nrows)
            break
        ## read next chunk
        vals <- tryCatch({
            read.delim(filecon, header=FALSE, sep="\t", skip=0, nrows=nrows,
                       stringsAsFactors=FALSE, quote="", colClasses=colClasses)
        }, error=function(err) {
            if (identical(conditionMessage(err), "no lines available in input"))
                data.frame()
            else stop(err)
        })
        overwrite <- FALSE
        append <- TRUE
    }
    close(filecon)
}

.indexTaxIds <-
    function(NCBIcon, tableName)
{
    ## create INDEX on tax_id
    sql <- paste0("
        CREATE INDEX IF NOT EXISTS ", tableName, "_taxId
        ON ", tableName, " (tax_id);")
    dbExecute(NCBIcon, sql)
}

.getFiles <-
    function(NCBIFilesDir, file, url, tableName)
{
    if (is.null(NCBIFilesDir)) {
        tmp <- tempfile()
        .tryDL(url, tmp)
    } else {
        tmp <- file.path(NCBIFilesDir, names(file))
        if (!file.exists(tmp))
            .tryDL(url, tmp)
    }
    tmp
}

## helpers for checking if NCBI 'cache database' is full (or not).
.setNCBIDateStamp  <-
    function(NCBIcon, tableName)
{
    tblNm <- paste0(tableName, '_date')
    vals = data.frame(date=as.character(Sys.Date()))
    dbWriteTable(NCBIcon, name=tblNm, value=vals, overwrite=TRUE)
}

.getNCBIDateStamp <-
    function(NCBIcon, tableName)
{
    tblNm <- paste0(tableName, '_date')
    dbGetQuery(NCBIcon, paste0("SELECT date FROM ",tblNm))$date
}

.isNCBICurrentWith <-
    function(NCBIcon, tableName)
{
    allTables <- dbListTables(NCBIcon)
    if (tableName %in% allTables) {
        DBdate <- as.character(.getNCBIDateStamp(NCBIcon, tableName))
        curdate <- as.character(Sys.Date())
        DBdate == curdate
    } else FALSE
}

## helper for checking if NCBI 'cache database' is populated yet
.isNCBIPopulatedWith <-
    function(NCBIcon, tableName)
{
    exists <- dbListTables(NCBIcon)
    tableName %in% exists
}

##########################################################
.downloadData <- function(file, tax_id, NCBIFilesDir, rebuildCache, verbose)
{
    if (verbose)
        message("getting data for ", names(file))

    ## NCBI connection
    if (is.null(NCBIFilesDir)) {
        NCBIcon <- dbConnect(SQLite(), dbname = tempfile())
        tmp <- tempfile()
    } else {
        NCBIcon <- dbConnect(SQLite(),
                             dbname = file.path(NCBIFilesDir, "NCBI.sqlite"))
        tmp <- file.path(NCBIFilesDir, names(file))
    }
    tableName <- sub(".gz","",names(file))

    ## download
    if (rebuildCache) {
        if(names(file) == "gene2unigene"){
            url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/", names(file))
        }else{
            url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/", names(file))
        }
        ## if DB table is not fresh OR if table is not populated
        if (!.isNCBICurrentWith(NCBIcon, tableName) ||
            !.isNCBIPopulatedWith(NCBIcon, tableName)) {
            if (verbose)
                message("rebuilding the cache")
            .tryDL(url, tmp)
            ## write to NCBI.sqlite db
            .writeToNCBIDB(NCBIcon, tableName, filepath=tmp, file)
            .setNCBIDateStamp(NCBIcon, tableName)
        }
    }
    

    if ("tax_id" %in% unlist(file)) {
        .indexTaxIds(NCBIcon, tableName)
    }
    

    ## get organism-specific data from NCBI.sqlite
    if (verbose)
        message("extracting data for our organism from : " ,tableName)
    if ("tax_id" %in% unlist(file)) { ## when there is a tax_id need to subset
        sql <- paste0("
            SELECT *
            FROM ", tableName, "
            WHERE ", tableName, ".tax_id = '", tax_id, "'")
        vals <- dbGetQuery(NCBIcon, sql)
        ## if things are empty: make sure we have correct number of cols
        if (dim(vals)[1] == 0) {## if after select there are no records here
            vals <- data.frame(t(seq_along(unlist(file))),
                               stringsAsFactors=FALSE)
            vals <- vals[FALSE,]
            colnames(vals) <- unlist(file)
        }
    } else { ## Just get the whole record
        message("getting all data for our organism from : " , tableName)
        sql <- paste0("SELECT * FROM ", tableName)
        vals <- dbGetQuery(NCBIcon, sql)
    }
    dbDisconnect(NCBIcon)
    ## remove row_names col
    vals[,!(colnames(vals) %in% 'row_names')]
}

## helper to populate base tables
.populateBaseTable <-
    function(con, sql, data, table)
{
    dbBegin(con)
    dbExecute(con, sql, unclass(unname(data)))
    dbCommit(con)
    message("table ", table, " filled")
}

## Convenience function for quickly making indices
.makeIndex <-
    function(con,table, field)
{
    indName <-paste("ind", table, field, sep="_")
    dbExecute(con, paste("
        CREATE INDEX IF NOT EXISTS", indName, "
        ON", table, "(", field, ")"))
}

## metadata stuff
.createMetadataTables <-
    function(con)
{
    dbExecute(con, "
        CREATE TABLE IF NOT EXISTS metadata (
            name VARCHAR(80) PRIMARY KEY,
            value VARCHAR(255)
        );")
    dbExecute(con, "
        CREATE TABLE IF NOT EXISTS map_metadata (
            map_name VARCHAR(80) NOT NULL,
            source_name VARCHAR(80) NOT NULL,
            source_url VARCHAR(255) NOT NULL,
            source_date VARCHAR(20) NOT NULL
         );")
    dbExecute(con, "
        CREATE TABLE IF NOT EXISTS map_counts (
            map_name VARCHAR(80) PRIMARY KEY,
            count INTEGER NOT NULL
        );")
}

.addMeta <-
    function(con, name, value)
{
    data = data.frame(name,value,stringsAsFactors=FALSE)
    sql <- "INSERT INTO metadata (name,value) VALUES(?,?)"
    .populateBaseTable(con, sql, data, "metadata")
}

.addMetadata <- function(con, tax_id, genus, species)
{
    name <- c("DBSCHEMAVERSION", "DBSCHEMA", "ORGANISM", "SPECIES", "CENTRALID",
              "Taxonomy ID",
              "EGSOURCEDATE", "EGSOURCENAME", "EGSOURCEURL",
              "GOSOURCEDATE", "GOSOURCENAME", "GOSOURCEURL",
              "GOEGSOURCEDATE", "GOEGSOURCENAME", "GOEGSOURCEURL",
              "Db type","Supporting package")
    value<- c("2.1","ORGANISM_DB", paste(genus, species), paste(genus, species),
              "EG",
              tax_id,
              date(), "Entrez Gene","ftp://ftp.ncbi.nlm.nih.gov/gene/DATA",
              .getGODate(),
              "Gene Ontology","ftp://ftp.geneontology.org/pub/go/godata",
              date(), "Entrez Gene","ftp://ftp.ncbi.nlm.nih.gov/gene/DATA",
              "OrgDb","AnnotationDbi")
    .addMeta(con, name, value)
}

.addMapMeta <-
    function(con, map_name, source_name, source_url, source_date)
{
    data = data.frame(map_name,source_name,source_url,
                      source_date,stringsAsFactors=FALSE)
    sql <- "
        INSERT INTO map_metadata (
            map_name, source_name, source_url, source_date
        ) VALUES(?,?,?,?);"
    .populateBaseTable(con, sql, data, "map_metadata")
}

.addMapMetadata <-
    function(con, tax_id, genus, species)
{
    map_name <- c("ENTREZID","GENENAME","SYMBOL","CHR","ACCNUM","REFSEQ","PMID",
                  "PMID2EG","UNIGENE","ALIAS2EG","GO2EG","GO2ALLEGS",
                  "GO","GO2ALLEGS")
    source_name <- c(rep("Entrez Gene",12),rep("Gene Ontology",2))
    source_url <- c(rep("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA",12),
                    rep("ftp://ftp.geneontology.org/pub/go/go",2))
    source_date <- c(rep(date(),12),rep(.getGODate(),2))
    .addMapMeta(con, map_name, source_name, source_url, source_date)

    ##clean up unwanted data
    ## if we used blast2GO, then we need to drop classic GO
    goSrcs <- dbGetQuery(con, "
        SELECT source_name
        FROM map_metadata
        WHERE map_name='GO2EG'")
    if ("Gene Ontology" %in% goSrcs & "blast2GO" %in% goSrcs) {
        dbExecute(con, "
            DELETE FROM map_metadata
            WHERE map_name = 'GO2EG'
                AND source_name='Gene Ontology'")
    }
}

## need to be able to generate the columns data from the files and the types.
.generateCols <-
    function (file)
{
    types <- c(tax_id = "INTEGER NOT NULL",
               gene_id = "INTEGER NOT NULL",
               go_id = "TEXT NOT NULL",
               evidence = "TEXT",
               go_qualifier = "TEXT",
               go_description = "TEXT",
               pubmed_id = "INTEGER",
               category = "TEXT",
               status = "TEXT",
               rna_accession = "TEXT",
               rna_gi = "INTEGER",
               protein_accession = "TEXT",
               protein_gi = "INTEGER",
               genomic_dna_accession = "TEXT",
               genomic_dna_gi = "INTEGER",
               genomic_start = "INTEGER",
               genomic_end = "INTEGER",
               orientation = "TEXT",
               assembly = "TEXT",
               unigene_id = "TEXT",
               symbol = "TEXT",
               locus_tag = "TEXT",
               synonyms = "TEXT",
               dbXrefs = "TEXT",
               chromosome = "TEXT",
               map_location = "TEXT",
               description = "TEXT",
               gene_type = "TEXT",
               nomenclature_symbol = "TEXT",
               nomenclature_name = "TEXT",
               nomenclature_status = "TEXT",
               other_designations = "TEXT",
               mim_id = "INTEGER NOT NULL",
               relation_type = "TEXT",
               refseq_id = "TEXT NOT NULL",
               uniprot_id = "TEXT NOT NULL")
    cols <- character()
    for(i in seq_len(length(file[[1]]))) {
        cols[i] <-  paste(file[[1]][i], types[file[[1]]][i])
    }
    paste(cols, collapse =",")
}

.mergeAndCleanAccession <-
    function(ac)
{
    ## pair 1st col with 2nd, clean and rename
    acr <- ac[ac[,2]!="-", 1:2]
    colnames(acr)[2] <- "accession"
    ## pair 1st col with 3rd, clean and rename
    acp <- ac[ac[,3]!="-",c(1,3)]
    colnames(acp)[2] <- "accession"
    ## merge these together
    ac <- rbind(acr,acp)
    ## Then clean off the period and trailing digits
    ac[,2] <- sub("\\.\\d+$","",ac[,2])
    ac
}

## also need somethign to generate the insert statement.
.generateTempTableINSERTStatement <-
    function (file)
{
    Qs <- paste(rep("?", length(file[[1]])),  collapse=",")
    paste0("
        INSERT INTO ", sub(".gz$", "", names(file)), "
        VALUES(", Qs, ");")
}


## need code to generate a table for each of the
## file below is like: files[1] or files[2]
.createTEMPNCBIBaseTable <-
    function (con, file, tax_id,
              NCBIFilesDir, rebuildCache, verbose)
{
    data <- .downloadData(file, tax_id, NCBIFilesDir=NCBIFilesDir,
                          rebuildCache=rebuildCache, verbose=verbose)
    table <- sub(".gz$","",names(file))
    if (is.null(table))
        stop("Unable to infer a table name.")
    cols <- .generateCols(file)
    sql<- paste("CREATE TABLE IF NOT EXISTS ", table, " (",cols,");")
    dbExecute(con, sql)

    message("Populating ", table," table:")
    sql <- .generateTempTableINSERTStatement(file)
    if (dim(data)[1] != 0) { ## ie. there has to be SOME data...
        .populateBaseTable(con, sql, data, table)
    } else {
        message("No data available for table ", table)
    }
}

## Download all the data files and make the cache DB
.setupBaseDBFromDLs <-
    function(files, tax_id, con, NCBIFilesDir, rebuildCache, verbose)
{
    for(i in seq_len(length(files))) {
        .createTEMPNCBIBaseTable(
            con, files[i], tax_id, NCBIFilesDir=NCBIFilesDir,
            rebuildCache=rebuildCache, verbose=verbose)
    }
    con
}

.dropOldTables <-
    function(con,fileNames)
{
    fileNames <- sub(".gz", "", fileNames)
    message("dropping table", fileNames)
    for(i in seq_len(length(fileNames))) {
        dbExecute(con, paste("DROP TABLE IF EXISTS", fileNames[i]))
    }
}


.generateOrgDbName <-
    function(genus, species)
{
    specStr <- paste0(toupper(substr(genus,1,1)), species)
    paste0("org.", specStr, ".eg")
}

.getGODate <-
    function()
{
    GOinf <- dbInfo(GO.db::GO_dbconn())
    GOinf[GOinf[["name"]]=="GOSOURCEDATE","value"]
}


## wanted when we want to know how many of something we can get
.computeSimpleEGMapCounts <-
    function(con, table, field)
{
    sql<- paste0("
        SELECT count(DISTINCT g.gene_id)
        FROM ", table, " AS t, genes AS g
        WHERE t._id=g._id AND t.", field, " NOT NULL")
    dbGetQuery(con,sql)[[1]]
}

## when we want to know how many of something there is (reverse mapping?)
.computeSimpleMapCounts <-
    function(con, table, field)
{
    sql<- paste0("
        SELECT COUNT(DISTINCT t.", field, ")
        FROM ", table, " AS t, genes AS g
        WHERE t._id=g._id AND t.", field, " NOT NULL")
    dbGetQuery(con,sql)[[1]]
}


## TODO: correct the counts as needed
.addMapCounts <-
    function(con, tax_id, genus, species)
{
    count <- c(
        GENENAME = .computeSimpleEGMapCounts(con, "gene_info", "gene_name"),
        SYMBOL = .computeSimpleEGMapCounts(con, "gene_info", "symbol"),
        SYMBOL2EG = .computeSimpleMapCounts(con, "gene_info", "symbol"),
        CHR = .computeSimpleEGMapCounts(con, "chromosomes", "chromosome"),
        ACCNUM = .computeSimpleEGMapCounts(con, "accessions", "accession"),
        REFSEQ = .computeSimpleEGMapCounts(con, "refseq", "accession"),
        REFSEQ2EG = .computeSimpleMapCounts(con, "refseq", "accession"),
        PMID = .computeSimpleEGMapCounts(con, "pubmed", "pubmed_id"),
        PMID2EG = .computeSimpleMapCounts(con, "pubmed", "pubmed_id"),
        UNIGENE = .computeSimpleEGMapCounts(con, "unigene", "unigene_id"),
        UNIGENE2EG = .computeSimpleMapCounts(con, "unigene", "unigene_id"),
        ALIAS2EG = .computeSimpleMapCounts(con, "alias", "alias_symbol"),
        GO = .computeSimpleEGMapCounts(con, "go", "go_id"),
        GO2EG = .computeSimpleMapCounts(con, "go", "go_id"),
        GO2ALLEGS = .computeSimpleMapCounts(con, "go_all", "go_id"),
        TOTAL = dbGetQuery(con,"SELECT count(DISTINCT gene_id) FROM genes")[[1]]
    )
    data = data.frame(
        map_names = names(count),
        count = unname(count),
        row.names = NULL
    )
    sql <- "INSERT INTO map_counts (map_name,count) VALUES(?,?)"
    .populateBaseTable(con, sql, data, "map_counts")
}

## list of primary NCBI files that we need to process
.primaryFiles <-
    function()
{
    list(gene2pubmed.gz=c("tax_id", "gene_id", "pubmed_id"),
         gene2accession.gz= c(
             "tax_id", "gene_id", "status", "rna_accession", "rna_gi",
             "protein_accession", "protein_gi", "genomic_dna_accession",
             "genomic_dna_gi", "genomic_start", "genomic_end", "orientation",
             "assembly", "peptide_accession", "peptide_gi", "symbol"),
         gene2refseq.gz=c(
             "tax_id", "gene_id", "status", "rna_accession", "rna_gi",
             "protein_accession", "protein_gi", "genomic_dna_accession",
             "genomic_dna_gi", "genomic_start", "genomic_end", "orientation",
             "assembly", "peptide_accession", "peptide_gi", "symbol"),
         gene2unigene=c("gene_id", "unigene_id"),
         gene_info.gz=c(
             "tax_id", "gene_id", "symbol", "locus_tag", "synonyms", "dbXrefs",
             "chromosome", "map_location", "description", "gene_type",
             "nomenclature_symbol", "nomenclature_name", "nomenclature_status",
             "other_designations", "modification_date", "feature_type"),
         gene2go.gz=c(
             "tax_id", "gene_id", "go_id", "evidence", "go_qualifier",
             "go_description", "pubmed_id", "category"))
}

#########################################################################
## Generate the database using the helper functions:
#########################################################################

makeOrgDbFromNCBI <-
    function(tax_id, genus, species, NCBIFilesDir, outputDir,
             rebuildCache=TRUE, verbose=TRUE)
{
    dbFileName <- paste0(.generateOrgDbName(genus,species), ".sqlite")
    dbFileName <- file.path(NCBIFilesDir, dbFileName)
    if (file.exists(dbFileName)) { file.remove(dbFileName) }
    con <- dbConnect(SQLite(), dbFileName)
    .createMetadataTables(con)  ## just makes the tables
    ## I need a list of files, along with their column names
    ## (needed for schema definitions later)
    ## IF ANY OF THESE gets moved in the source files
    ## then things will be in the wrong place!
    files = .primaryFiles()
    .setupBaseDBFromDLs(files, tax_id, con, NCBIFilesDir=NCBIFilesDir,
                        rebuildCache=rebuildCache, verbose=verbose)

    ## Add metadata:
    .addMetadata(con, tax_id, genus, species)
    ## Add map_metadata:
    .addMapMetadata(con, tax_id, genus, species)

    ## Make the central table
    egs <- dbGetQuery(con, "
        SELECT distinct gene_id FROM gene_info")[,1]
    .makeCentralTable(egs, con)

    ## Make the other tables:
    ## gene_info
    symbs <- dbGetQuery(con, "
        SELECT distinct gene_id, description, symbol FROM gene_info")
    colnames(symbs) <- c("gene_id", "gene_name", "symbol")
    ## constraint requires that we always have BOTH a symbol and a name.
    symbs <- symbs[!is.na(symbs[,2]) & !is.na(symbs[,3]),] # strict constraint!
    ## Only 10 things fail this strict constraint for human.
    ## ALL because of no gene symbol (but that have a name)
    ## TODO: fix this so that I can make this table without
    ## such a strict constraint
    .makeSimpleTable(symbs, table="gene_info_temp", con,
                     fieldNameLens=c(255,80))

    ## alias ## requires sub-parsing.
    alias <- dbGetQuery(con, "
        SELECT distinct gene_id, synonyms FROM gene_info")
    aliases <- sapply(alias[,2], strsplit, "\\|")
    numAlias <- sapply(aliases, length)
    alGenes <- rep(alias[,1], numAlias)
    alias <- data.frame(gene_id=alGenes,alias_symbol=unlist(aliases))
    .makeSimpleTable(alias, table="alias", con)

    ## chr
    chrs <- dbGetQuery(con, "
        SELECT distinct gene_id, chromosome FROM gene_info")
    .makeSimpleTable(chrs, table="chromosomes", con)


    ## pubmed
    pm <- dbGetQuery(con, "
        SELECT distinct gene_id, pubmed_id FROM gene2pubmed")
    .makeSimpleTable(pm, table="pubmed", con)

    ## refseq ## requires sub-parsing.
    rs <- dbGetQuery(con, "
        SELECT DISTINCT gene_id, rna_accession, protein_accession
        FROM gene2refseq")
    rs <- .mergeAndCleanAccession(rs)
    .makeSimpleTable(rs, table="refseq", con)

    ## accessions
    ac <- dbGetQuery(con, "
        SELECT DISTINCT gene_id, rna_accession, protein_accession
        FROM gene2accession")
    ac <- .mergeAndCleanAccession(ac)
    .makeSimpleTable(ac, table="accessions", con)

    ## unigene
    ug <- dbGetQuery(con, "
        SELECT distinct g.gene_id, u.unigene_id
        FROM gene2unigene AS u, gene_info AS g
        WHERE u.gene_id=g.gene_id")
    .makeSimpleTable(ug, table="unigene", con)

    ## Make the GO tables:
    .makeGOTablesFromNCBI(con)

    ## Drop all the older tables (which will include the original "gene_info").
    .dropOldTables(con,names(files))
    ## Rename "gene_info_temp" to be just "gene_info":
    dbExecute(con,"ALTER TABLE gene_info_temp RENAME TO gene_info")
    ## add GO views to the DB
    makeGOViews(con)

    ## Add map_counts: (must happen at end)
    .addMapCounts(con, tax_id, genus, species)
    dbDisconnect(con)

}

## For argument checking:
.isSingleString <-
    function(x)
{
    is.atomic(x) && length(x) == 1L && is.character(x)
}
.isSingleStringOrNA <-
    function(x)
{
    is.atomic(x) && length(x) == 1L && (is.character(x) || is.na(x))
}
.isSingleStringOrNull <-
    function(x)
{
    (is.atomic(x) && length(x) == 1L && is.character(x)) ||
    (is.atomic(x) && length(x) == 0L && is.null(x))
}

## OLDER function to make the package:
OLD_makeOrgPackageFromNCBI <-
    function(version, maintainer, author, outputDir, tax_id, genus, species,
             NCBIFilesDir, rebuildCache)
{
    message(
        "If this is the 1st time you have run this function, it may take, a ",
        "long time (over an hour) to download needed files and assemble a ",
        "33 GB cache databse in the NCBIFilesDir directory.  Subsequent calls ",
        "to this function should be faster (seconds).  The cache will try to ",
        "rebuild once per day.")

    ## Arguement checking:
    if (!.isSingleString(version))
        stop("'version' must be a single string")
    if (!.isSingleString(maintainer))
        stop("'maintainer' must be a single string")
    if (!.isSingleString(author))
        stop("'author' must be a single string")
    if (outputDir!="." && file.access(outputDir)[[1]]!=0) {
        stop("Selected outputDir '", outputDir,"' does not exist.")}
    if(!(isSingleNumber(tax_id) || .isSingleString(tax_id)))
        stop("'tax_id' must be a single integer")
    if (!is.integer(tax_id))
        tax_id <- as.integer(tax_id)
    if (!.isSingleString(genus))
        stop("'genus' must be a single string")
    if (!.isSingleString(species))
        stop("'species' must be a single string")
    if (!.isSingleStringOrNull(NCBIFilesDir))
        stop("'NCBIFilesDir' argument needs to be a single string or NULL")

    makeOrgDbFromNCBI(tax_id=tax_id, genus=genus, species=species,
                      NCBIFilesDir=NCBIFilesDir, outputDir, rebuildCache)

    dbName <- .generateOrgDbName(genus,species)
    dbfile <- paste0(dbName, ".sqlite")

    seed <- new("AnnDbPkgSeed",
                Package= paste0(dbName, ".db"),
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
    ## return handle to the DB file name
    dbfile
}



#############################################################################
## Code specific for the new process.

## based on tax_id do we need to look at unigene?
.doWeNeedUnigene <-
    function(tax_id)
{
    file <- system.file('extdata','unigene.txt',package='AnnotationForge')
    res <- read.table(file, skip=1, sep=":")
    tax_id %in% res[[2]]
}

.makeBaseDBFromDLs <-
    function(files, tax_id, con, NCBIFilesDir, rebuildCache, verbose)
{
    ## Test if we need to bother with unigene
    if (!.doWeNeedUnigene(tax_id))
        files <- files[!(names(files) %in% 'gene2unigene')]
    ## get/update the data
    data <- list(length(files))
    if (verbose)
        message("starting download for ", sprintf("\n[%d] %s",
                                                  seq_along(files),
                                                  names(files)))
    for(i in seq_len(length(files))) {
        res <- .downloadData(files[i], tax_id, NCBIFilesDir=NCBIFilesDir,
                             rebuildCache=rebuildCache, verbose=verbose)
        data[[i]] <- res
    }
    names(data) <- names(files)
    data
}

## Helper functions for extracting data from a frame, making sure
## there *is* some data and renaming it etc.
.tossHyphens <-
    function(df)
{
    cols <- seq_len(ncol(df)) ## get range of columns in df
    cols <- cols[-1] ## drop the 1st one
    for(i in cols) {
        df <- df[df[[i]]!='-',]
    }
    df
}

## extraction function
.extract <-
    function(rawData, keepTable, keepCols)
{
    rawData[[keepTable]][,keepCols]
}

## rename and alter the data list (if needed)
.rename <-
    function(data, dataSub, newName, newCols)
{
    dataSub <- .tossHyphens(dataSub)
    if (dim(dataSub)[1] > 0) { ## if there are rows...
        data[[newName]] <- unique(dataSub)
        colnames(data[[newName]]) <- newCols
    }
    data
}

## extract AND rename etc.
.procDat <-
    function(rawData, keepTable, keepCols, data, newName, newCols)
{
    pubSub <- .extract(rawData, keepTable, keepCols)
    .rename(data, pubSub, newName, newCols)
}

## newer .getBlast2GO is based (very loosely) on older .Blast2GOData()
## Alternative Blast2GO (refactor in progress)
## The file we have to proces has specs here:
## http://www.b2gfar.org/fileformat

## need a version of this that does not throw in towel if GO is missing
.tryDL2 <-
    function(url, tmp, times = 2)
{
    for(i in seq_len(times)) {
        ## tryResult <- try( download.file(url, tmp, quiet=TRUE) , silent=TRUE)
        tryResult <- try( .downloadAndSaveToTemp(url, tmp) , silent=TRUE)
        if (is(tryResult,"try-error") && i < times) {
            Sys.sleep(20)
        } else if (is(tryResult,"try-error") && i >= times) {
            msg = c("url access failed after ", times, " attempts; url: ", url)
            warning(paste(strwrap(msg,exdent=2), collapse="\n"))
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
}

.emptyGoFrame <-
    function()
{
    data.frame(gene_id=character(), go_id=character(),
               evidence=character(), stringsAsFactors=FALSE)
}

.getBlast2GO <-
    function(tax_id, refseq, accs)
{
    url = paste0(
        "http://www.b2gfar.org/_media/species:data:", tax_id,".annot.zip")
    tmp <- tempfile()
    ##if we can download the file, process it.
    if (.tryDL2(url,tmp)) {
        rawVals <- read.delim(unzip(tmp), header=FALSE, sep="\t", quote="",
                              stringsAsFactors=FALSE)
    } else {## return empty thing right now.
        return(.emptyGoFrame())
    }
    ## if there are no refseqs and now accessions either (from NCBI) -
    ## then bail
    if (is.null(refseq) && is.null(accs)) {
        warning("No GO data was able to be matched to this tax ID.")
        return(.emptyGoFrame())
    }
    ## I will need to so extra stuff here to match up categories etc.
    ## (vals has to look like gene2go would, and I have to join to
    ## refseq and to accession just to get my EGs)
    RSVals <- rawVals[grep("refseq",rawVals[,4]),c(2,3)]
    GBVals <- rawVals[grep("genbank",rawVals[,4]),c(2,6)]
    colnames(GBVals) <- colnames(RSVals) <- c("go_id","accession")

    ## Get gene_ids matches to refseq and merge
    RSIDs <- merge(RSVals, refseq, by.x='accession', by.y='REFSEQ')
    ## Get gene_ids matches to accessions and merge
    GBIDs <- merge(GBVals, accs, by.x='accession', by.y='ACCNUM')

    ## Combine refseq matches and accession mathes
    vals <- unique(rbind(GBIDs, RSIDs))
    rownames(vals) <- NULL

    ## Then do this to make the final data.frame:
    data.frame(gene_id = as.character(vals[["GID"]]),
               go_id = vals[["go_id"]],
               evidence = rep("IEA", length(vals[["go_id"]])),
               stringsAsFactors=FALSE)
}

## Helper to get alternate GO data from UniProt and add it to the NCBI DB.
.downloadAndPopulateAltGOData <-
    function(NCBIcon, NCBIFilesDir, rebuildCache)
{
    dest <- file.path(NCBIFilesDir, "idmapping_selected.tab.gz")
    if (rebuildCache) {
        #  This url has been flaky in the past
        #  See https://www.uniprot.org/downloads
        #  Troublshooting in the past involved temporarily changing this url
        #     to use the https protcol url:
        #     https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
        #url <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
        url <- "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
        ## RCurl has TLS problems on Windows use download.file instead
        ## loadNamespace("RCurl")
        ## f <- RCurl::CFILE(dest, mode="wb")
        ## RCurl::curlPerform(url = url, writedata = f@ref)
        options(timeout = max(1000, getOption("timeout")))
        status <- download.file(url, dest, quiet = TRUE)
        if(status != 0L) stop("Download of alternative GO data failed. Please retry!\n", call. = FALSE)
    }
    ## create table and set up indices
    dbExecute(NCBIcon, "
        CREATE TABLE IF NOT EXISTS altGO (
            EntrezGene TEXT NOT NULL,
            GO TEXT NOT NULL,
            NCBItaxon TEXT NOT NULL
        )")
    dbExecute(NCBIcon, "CREATE INDEX IF NOT EXISTS geneidAltGO on altGO(EntrezGene)")
    dbExecute(NCBIcon, "CREATE INDEX IF NOT EXISTS taxidAltGO on altGO(NCBItaxon)")
    ## Select fields 3 (EntrezGene), 7(GO) and 13 (NCBI-taxon)
    colClasses <- c("NULL","NULL","character","NULL","NULL","NULL",
                    "character","NULL","NULL","NULL","NULL","NULL",
                    "character","NULL","NULL","NULL","NULL","NULL",
                    "NULL","NULL","NULL","NULL")
    sql <- "INSERT INTO altGO (EntrezGene,GO,NCBItaxon) VALUES(?,?,?)"

    ## write to table
    .writeToNCBIDB2(NCBIcon, sql=sql, filePath=dest, tableName="altGO",
                    colNames=c("EntrezGene", "GO", "NCBItaxon"),
                    colClasses=colClasses)
    ##data <-  read.delim(dest, header=FALSE, sep="\t", quote="",
    ##                    stringsAsFactors=FALSE, colClasses=colClasses)
    ##.populateBaseTable(NCBIcon, sql, data, "metadata")
    ## And don't forget to stamp out a date table
    .setNCBIDateStamp(NCBIcon, 'altGO')
}

## FIXME: merge with .writeToNCBIDB
.writeToNCBIDB2 <-
    function(NCBIcon, sql, filePath, tableName, colNames, colClasses,
             initialSkip=0)
{
    nrows <- 1000000
    overwrite <- TRUE
    append <- FALSE
    ## first chunk
    filecon <- file(description=filePath, open="r")
    vals <- read.delim(filecon, header=FALSE, sep="\t", skip=initialSkip,
                       nrows=nrows,
                       stringsAsFactors=FALSE, quote="", colClasses=colClasses)
    repeat {
        if (nrow(vals) == 0)
            break
        ## process
        if ('tax_id' %in% colnames(vals))
            vals[['tax_id']] <- as.integer(vals[['tax_id']])
        colnames(vals) <- colNames
        dbWriteTable(NCBIcon, tableName, vals,
                     overwrite=overwrite, append=append)
        ## break if last chunk was final chunk
        if (nrow(vals) != nrows)
            break
        ## read next chunk
        vals <- tryCatch({
            read.delim(filecon, header=FALSE, sep="\t", skip=0, nrows=nrows,
                       stringsAsFactors=FALSE, quote="", colClasses=colClasses)
        }, error=function(err) {
           if (identical(conditionMessage(err), "no lines available in input"))
              data.frame()
           else stop(err)
        })
        overwrite <- FALSE
        append <- TRUE
    }
    close(filecon)
}

## Helper to allow splitting of a two collumn data frame where
## either column may have "; " separated values
splitBy <-
    function(data, splitCol=1)
{
    splits <- strsplit(data[[splitCol]],split="; ")
    splitLens <- unlist(lapply(splits, length))

    if (splitCol==1) {
        dups <- rep(data[[2]], times=splitLens)
        data <- data.frame(eg_id=unlist(splits), go_id=dups,
                           stringsAsFactors=FALSE)
    }
    if (splitCol==2) {
        dups <- rep(data[[1]], times=splitLens)
        data <- data.frame(eg_id=dups, go_id=unlist(splits),
                           stringsAsFactors=FALSE)
    }
    data
}

## GO Data for organisms with no data at NCBI
## FIXME: set a timeout? (2.8GB file)
.getAltGOData <-
    function(NCBIcon, NCBIFilesDir, tax_id, rebuildCache)
{
    ## First get the data and populate it (if necessary)
    if (!.isNCBICurrentWith(NCBIcon, 'altGO') ||
       !.isNCBIPopulatedWith(NCBIcon, 'altGO')) {
        .downloadAndPopulateAltGOData(NCBIcon, NCBIFilesDir, rebuildCache)
    }
    
    ## Then get out the data that we *want* using another query with the tax_id
    sql <- paste0("
        SELECT EntrezGene, GO
        FROM altGO
        WHERE NCBItaxon = '", tax_id, "'")
    res <- dbGetQuery(NCBIcon, sql)
    ## THEN, we have to unpack the GO terms (expand them)
    ## GOs <- strsplit(res$GO,split="; ")
    ## goLens <- unlist(lapply(GOs, length))
    ## entrez <- rep(res$EntrezGene, times=goLens)
    ## evidence <- rep("IEA", times=length(entrez))
    ## data.frame(gene_id=entrez, go_id=unlist(GOs), evidence=evidence,
    ##            stringsAsFactors=FALSE)

    ## DROP records with no entrez gene ID or no GO ID
    res <- res[!(res$EntrezGene == "" | res$GO == ""), ]
    if (nrow(res) == 0)
        return(data.frame(gene_id="", go_id="", evidence="",
                          stringsAsFactors=FALSE)[FALSE,])

    ## split two column data frames at will...
    data <- splitBy(res, splitCol=2)
    ## then split it the other way too fully expand it:
    data <- splitBy(data, splitCol=1)
    ## enforce stringsAsFactors=FALSE
    data.frame(gene_id=data[[1]], go_id=data[[2]], evidence="IEA",
               stringsAsFactors=FALSE)
}


## Helper to get complete list of entrez gene IDs
.mkGIdQ <-
    function(tax_id, tableName)
{
    paste0("
        SELECT gene_id
        FROM ", tableName, "
        WHERE ", tableName, ".tax_id = '", tax_id, "'")
}
.asV <-
    function(df)
{
    unique(as.character(t(df)))
}
.getAllEnrezGeneIdsFromNCBI <-
    function(tax_id, NCBIcon)
{
    id1 <- .asV(dbGetQuery(NCBIcon, .mkGIdQ(tax_id, 'gene_info')))
    id2 <- .asV(dbGetQuery(NCBIcon, .mkGIdQ(tax_id, 'gene2refseq')))
    id3 <- .asV(dbGetQuery(NCBIcon, .mkGIdQ(tax_id, 'gene2go')))
    id4 <- .asV(dbGetQuery(NCBIcon, .mkGIdQ(tax_id, 'gene2pubmed')))
    id5 <- .asV(dbGetQuery(NCBIcon, .mkGIdQ(tax_id, 'gene2accession')))
    unique(c(id1,id2,id3,id4,id5))
}

## Create list of data frames from the NCBI data
prepareDataFromNCBI <-
    function(tax_id, NCBIFilesDir, outputDir, rebuildCache, verbose, ensemblVersion)
{
    ## Files that need to be processed
    files <- .primaryFiles()
    NCBIcon <- dbConnect(SQLite(),
                         dbname = file.path(NCBIFilesDir, "NCBI.sqlite"))
    rawData <- .makeBaseDBFromDLs(files, tax_id, NCBIcon, NCBIFilesDir,
                                  rebuildCache, verbose)
    data <- list()
    ## pubmed
    if (verbose)
        message("processing gene2pubmed")
    data <- .procDat(rawData,
                     keepTable='gene2pubmed.gz',
                     keepCols=c('gene_id','pubmed_id'),
                     data,
                     newName='pubmed',
                     newCols=c('GID', 'PMID'))
    ## chromosomes
    if (verbose)
        message("processing gene_info: chromosomes")
    data <- .procDat(rawData,
                     keepTable='gene_info.gz',
                     keepCols=c('gene_id','chromosome'),
                     data,
                     newName='chromosomes',
                     newCols=c('GID', 'CHR'))
    ## gene_info
    if (verbose)
        message("processing gene_info: description")
    data <- .procDat(rawData,
                     keepTable='gene_info.gz',
                     keepCols=c('gene_id','description','symbol'),
                     data,
                     newName='gene_info',
                     newCols=c('GID', 'GENENAME','SYMBOL'))

    if (!length(data))
        stop("no information found for species with tax id ", tax_id)

    ## entrez genes
    ## FIXME: look at this
    egIDs <- .getAllEnrezGeneIdsFromNCBI(tax_id, NCBIcon)
    egData <- data.frame(GID=egIDs, ENTREZID=egIDs,stringsAsFactors=FALSE)
    data <- .rename(data,
                    egData,
                    newName='entrez_genes',
                    newCols=c('GID', 'ENTREZID'))

    ## Alias table:
    ## Must contain both symbols AND info from synonyms.
    ## For alias I need to massage the synonyms field from gene_info
    ## and add it to symbols from the same.
    if (verbose)
        message("processing alias data")
    aliasSQL <- paste0("
        SELECT distinct gene_id, synonyms
        FROM gene_info
        WHERE tax_id ='", tax_id, "'")
    alias <- dbGetQuery(NCBIcon,aliasSQL)
    aliases <- sapply(alias[,2], strsplit, "\\|")
    numAlias <- sapply(aliases, length)
    allGenes <- rep(alias[,1], numAlias)
    aliasDat <- data.frame(gene_id=allGenes,alias_symbol=unlist(aliases),
                           stringsAsFactors=FALSE)
    symbolDat <- .extract(rawData,
                          keepTable='gene_info.gz',
                          keepCols=c('gene_id','symbol'))
    colnames(aliasDat) <- NULL
    colnames(symbolDat) <- NULL
    if (nrow(aliasDat) == 0L)
        faliasDat <- as.matrix(symbolDat)
    else
        faliasDat <- unique(rbind(as.matrix(aliasDat), as.matrix(symbolDat)))
    rownames(faliasDat) <- NULL
    data <- .rename(data,
                    as.data.frame(faliasDat, stringsAsFactors=FALSE),
                    newName='alias',
                    newCols=c('GID', 'ALIAS'))

    ## refseq:
    ## Custom handling to combine a couple of columns into one.
    if (verbose)
        message("processing refseq data")
    refDat1 <- .extract(rawData,
                     keepTable='gene2refseq.gz',
                     keepCols=c('gene_id','protein_accession'))
    refDat2 <- .extract(rawData,
                     keepTable='gene2refseq.gz',
                     keepCols=c('gene_id','rna_accession'))
    colnames(refDat1) <- NULL
    colnames(refDat2) <- NULL
    refDat <- rbind(as.matrix(refDat1),as.matrix(refDat2))
    rownames(refDat) <- NULL
    data <- .rename(data,
                    as.data.frame(refDat, stringsAsFactors=FALSE),
                    newName='refseq',
                    newCols=c('GID', 'REFSEQ'))

    ## accession:
    ## Similar to refseq
    if (verbose)
        message("processing accession data")
    accDat1 <- .extract(rawData,
                     keepTable='gene2accession.gz',
                     keepCols=c('gene_id','protein_accession'))
    accDat2 <- .extract(rawData,
                     keepTable='gene2accession.gz',
                     keepCols=c('gene_id','rna_accession'))
    colnames(accDat1) <- NULL
    colnames(accDat2) <- NULL
    accDat <- rbind(as.matrix(accDat1),as.matrix(accDat2))
    rownames(accDat) <- NULL
    data <- .rename(data,
                    as.data.frame(accDat, stringsAsFactors=FALSE),
                    newName='accessions',
                    newCols=c('GID', 'ACCNUM'))

    ## GO:
    if (verbose)
        message("processing GO data")
    ## Step I: try gene2go
    goDat <- .extract(rawData,
                      keepTable='gene2go.gz',
                      keepCols=c('gene_id','go_id','evidence'))
    ## Step II: try Blast2GO
    if (dim(goDat)[1] == 0) {
        goDat <- .getAltGOData(NCBIcon, NCBIFilesDir, tax_id, rebuildCache)
    }
    ## .rename checks to make sure there actually ARE GO IDs
    data <- .rename(data,
                    goDat,
                    newName='go',
                    newCols=c("GID","GO","EVIDENCE"))

    ## If we have ensembl gene ids as data, then lets include it too
    if (tax_id %in% names(available.ensembl.datasets(tax_id, release=ensemblVersion))) {
        if (verbose)
            message("processing ensembl gene id data")
        ensDat <- .getEnsemblData(taxId=tax_id, release=ensemblVersion)
        data <- .rename(data,
                        ensDat,
                        newName='ensembl',
                        newCols=c("GID","ENSEMBL"))
    }

    ## unigene:
    ## Custom extraction of relevant genes ONLY (do this last)
    if (.doWeNeedUnigene(tax_id)) {
        if (verbose)
            message("processing unigene data")
        uniDat <- .extract(rawData,
                           keepTable='gene2unigene',
                           keepCols=c('gene_id','unigene_id'))
        ## Filter things not in our set of entrez gene IDs
        uniDat <- uniDat[uniDat$gene_id %in% egIDs,]
        data <- .rename(data,
                        dataSub=uniDat,
                        newName='unigene',
                        newCols=c('GID', 'UNIGENE'))
    }
    dbDisconnect(NCBIcon)
    data
}

## For non-traditional packages, i.e., those that appear only in AnnotationHub
NEW_makeOrgPackageFromNCBI <-
    function(version, maintainer, author, outputDir, tax_id, genus, species,
             NCBIFilesDir, databaseOnly, rebuildCache, verbose, ensemblVersion)
{
    if (rebuildCache)
        message("If files are not cached locally this may take ",
                "awhile to assemble a 33 GB cache databse in the ",
                "NCBIFilesDir directory. Subsequent calls to this ",
                "function should be faster (seconds). The cache will ",
                "try to rebuild once per day.",
                "Please also see AnnotationHub for some pre-built",
                "OrgDb downloads")
    
    

    if (!.isSingleString(version))
        stop("'version' must be a single string")
    if (!.isSingleString(maintainer))
        stop("'maintainer' must be a single string")
    if (!.isSingleString(author))
        stop("'author' must be a single string")
    if (outputDir!="." && file.access(outputDir)[[1]]!=0) {
        stop("Selected outputDir '", outputDir,"' does not exist.")}
    if(!(isSingleNumber(tax_id) || .isSingleString(tax_id)))
        stop("'tax_id' must be a single integer")
    if (!is.integer(tax_id))
        tax_id <- as.integer(tax_id)
    if (!.isSingleStringOrNull(genus))
        stop("'genus' must be a single string or NULL")
    if (!.isSingleStringOrNull(species))
        stop("'species' must be a single string or NULL")
    if (!.isSingleStringOrNull(NCBIFilesDir))
        stop("'NCBIFilesDir' argument needs to be a single string or NULL")

    ## genus and species
    if (is.null(genus))
        genus <- GenomeInfoDb:::lookup_organism_by_tax_id(tax_id)[['genus']]
    if (is.null(species)) {
        species <- GenomeInfoDb:::lookup_organism_by_tax_id(tax_id)[['species']]
        species <- gsub(' ', '.', species)
    }

    if (verbose)
        message("preparing data from NCBI ...")
    data <- prepareDataFromNCBI(tax_id, NCBIFilesDir, outputDir,
                                rebuildCache, verbose, ensemblVersion=ensemblVersion)
    dbName <- .generateOrgDbName(genus,species)
    dbfile <- paste0(dbName, ".sqlite")

    ## If there is go data use the goTable argument.
    goTable <- NA
    if (!is.null(data[['go']]))
        goTable <- "go"

    if (verbose)
        message("making the OrgDb package ...")
    .makeOrgPackage(data,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    outputDir=outputDir,
                    tax_id=tax_id,
                    genus=genus,
                    species=species,
                    goTable=goTable,
                    databaseOnly=databaseOnly, verbose=verbose)

    ## return handle to the DB name
    if (verbose)
        message("complete!")
    dbfile
}

## function wrapper to make the package:
makeOrgPackageFromNCBI <-
    function(version, maintainer, author, outputDir=getwd(), tax_id,
             genus=NULL, species=NULL, NCBIFilesDir=getwd(), databaseOnly=FALSE,
             useDeprecatedStyle=FALSE, rebuildCache=TRUE, verbose=TRUE, ensemblVersion=NULL)
{
    if (useDeprecatedStyle==TRUE) {
        dbname <- OLD_makeOrgPackageFromNCBI(version,maintainer,author,
                                             outputDir,tax_id,genus,
                                             species,NCBIFilesDir,
                                             rebuildCache=rebuildCache)
    } else {
        dbname <- NEW_makeOrgPackageFromNCBI(version,maintainer,author,
                                             outputDir,tax_id,genus,species,
                                             NCBIFilesDir,databaseOnly,
                                             rebuildCache=rebuildCache,
                                             verbose=verbose,
                                             ensemblVersion=ensemblVersion)
    }
    ## return handle to the db name
    dbname
}

################################################################################
################################################################################
## Code for adding Ensembl IDs for those species supported by ensembl
## with GTF files.

## use the Ensembl Rest API to get the current Ensembl version
getCurrentEnsemblRelease <- function() {
    
    res <- httr::GET("https://rest.ensembl.org/info/data/", httr::accept_json())
    version <- as.integer(unlist(httr::content(res)))
    
    return(version)
}

## list files in Ensembl FTP TSV directory
## this is used to identify if there is an ensembl<->entrez mapping file
listFilesInEnsemblFTP <- function(species, release, dir) {
    
    loadNamespace("curl")
    
    ftp_dir <- sprintf("ftp://ftp.ensembl.org/pub/release-%s/%s/%s/", release, dir, species)
    
    list_files <- curl::new_handle()
    curl::handle_setopt(list_files, ftp_use_epsv = TRUE, dirlistonly = TRUE)
    con <- curl::curl(url = ftp_dir, open = "r", handle = list_files)
    files <- readLines(con)
    close(con)
    
    return(files)
    
}

## lets make a helper to troll the ftp site and get ensembl to entrez
## gene ID data.  The names will be tax ids...
getFastaSpeciesDirs <-
    function(release=NULL)
    {
        if (is.null(release)) {
            release <- getCurrentEnsemblRelease()
        }

        listing <- listFilesInEnsemblFTP(release = release, species = "", dir = "mysql")
        res <- grep(pattern = paste0("_core_", release, "_"), listing, value = TRUE)
        
        specNames <- available.FastaEnsemblSpecies(res, release=release)
        taxdb <- GenomeInfoDb::loadTaxonomyDb()
        organisms <- trimws(paste(taxdb[["genus"]], taxdb[["species"]]))
        idx <- match(specNames, organisms)
        res <- res[!is.na(idx)]
        names(res) <- as.integer(taxdb[["tax_id"]][idx[!is.na(idx)]])
        res
    }


## Helper for getting precise genus and species available via
## FTP. (and their Tax IDs)
available.FastaEnsemblSpecies <-
    function(speciesDirs=NULL, release=NULL)
{
    if (is.null(speciesDirs)) {
        speciesDirs <- getFastaSpeciesDirs(release=release)
    }
    species <- sub('_core_\\d+_\\d+','',speciesDirs)
    genSpec <- gsub('_', ' ', species)
    .upperCaseGenus <- function(str) {
        paste0(toupper(substr(str,start=1,stop=1)),
               substr(str,start=2,stop=nchar(str)))
    }
    unlist(lapply(genSpec, .upperCaseGenus))
}
## Basically this is so that I can look up tax IDs like this:
## library(GenomeInfoDb); specs <- available.FastaEnsemblSpecies();
## lapply(specs, GenomeInfoDb:::lookup_tax_id_by_organism)

## helper for parsing strings into g.species format
g.species <-
    function(str)
{
    if (str %in% c("Canis familiaris", "Canis lupus familiaris"))
        return("clfamiliaris")
    strVec <- unlist(strsplit(str, split=' '))
    firstLetter <- tolower(substr(strVec[1],start=1,stop=1))
    theLast <- strVec[length(strVec)]
    paste0(firstLetter, theLast)
}

## helper to make sure that we *have* entrez gene IDs to map to at ensembl!
.ensemblMapsToEntrezId <-
    function(taxId, species, release=NULL)
{
    message("TaxID: ", taxId)

    files <- listFilesInEnsemblFTP(species = species, release = release, dir = 'tsv')
    
    any(grepl('entrez.tsv.gz$', files))
}

## the available.ensembl.datasets function takes 20 seconds to make a small
## vector.  So stash the results here for faster reuse/access on
## subsequent calls
## ## TODO: this caching mechanism may not be necessary anymore with the move away from biomaRt
## ## Nonetheless it still works so doesn't need to be removed - MLSmith 2024-04-08
ensemblDatasets <- new.env(hash=TRUE, parent=emptyenv())
## populate the above with: available.ensembl.datasets()

## I want to get the availble datasets for all the requested fasta species.
available.ensembl.datasets <-
    function(taxId=NULL, release=NULL)
{
  ## if the package wide vector is not empty, return it now.
    if (exists('ensDatSets', envir=ensemblDatasets)) {
        datSets <- get('ensDatSets', envir=ensemblDatasets)
        if (!is.null(taxId)) {
            taxId <- setdiff(taxId, names(datSets))
        }
    } else {
        datSets <- NULL
    }
    if (is.null(datSets) || is.null(taxId) || length(taxId)>0) {
        
        if (is.null(release)) {
            release <- getCurrentEnsemblRelease()
        }
        
        fastaSpecs <- available.FastaEnsemblSpecies(release=release)
        fasta_specs <- tolower(gsub(x = fastaSpecs, pattern = " ", replacement = "_"))
        if (is.null(taxId)) {
            taxId <- setdiff(names(fastaSpecs), names(datSets))
        }
        g.specs <- unlist(lapply(fastaSpecs, g.species))
        names(fasta_specs) <- names(g.specs)
        neededDatSets <- fasta_specs[names(fasta_specs) %in% taxId]
        if (length(neededDatSets)>0) {
            if (length(neededDatSets)>4) {
                ## Friendly message
                message(wmsg("Please be patient while we work out which organisms can
                              be annotated with ensembl IDs. "))
            }
            ## Remove dataSets that don't map to EntrezIds:
            legitTaxIdxs <- mapply(.ensemblMapsToEntrezId, 
                                   names(neededDatSets), neededDatSets, 
                                   MoreArgs = list(release = release))
            neededDatSets <- neededDatSets[legitTaxIdxs]
            datSets <- c(datSets, neededDatSets)

            assign("ensDatSets", datSets,
                   envir=ensemblDatasets)
        }
    }
    datSets
}

## ## now get those from the ensembl marts
.getEnsemblData <-
    function(taxId, release=NULL)
{
    if (is.null(release)) {
        release <- getCurrentEnsemblRelease()
    }
        
    datSets <- available.ensembl.datasets(taxId, release=release)
    datSet <- datSets[names(datSets) %in% taxId]
    if(length(datSet) == 0) {
        ## TODO: Is this a sensible return value?  Should it be an error?
        stop('No Ensembl dataset found matching the taxId: ', taxId)
    } else {
        ## list the files in the species level TSV folder.  There should always be a entrez mapping
        ## file as we checked for it earlier with available.ensembl.datasets()
        ## Find it and construct a URL with that file name.
        files <- listFilesInEnsemblFTP(species = datSet, release = release, dir = 'tsv')
        entrez_map_file <- grep(files, pattern = 'entrez.tsv.gz$', value = TRUE)
        url <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/tsv/%s/%s', release, datSet, entrez_map_file)
        
        ## download the mapping file, extract the relevant columns, rename, and return unique rows
        dest_file <- tempfile()
        download.file(url, destfile = dest_file, mode = "wb")
        res <- read.delim(dest_file)[,c('xref', 'gene_stable_id')]
        colnames(res) <- c("gene_id","ensembl")
        res <- res[!is.na(res$gene_id),]
        unique(res)
    }

}
