#########################################################################
#########################################################################
####         Helper functions for generating a generic DB            ####
#########################################################################
#########################################################################

## This makes the central table of an EG DB.
.makeCentralTable <- function(entrez, con){
  message("Populating genes table:")
  sql<- paste("    CREATE TABLE IF NOT EXISTS genes (
      _id INTEGER PRIMARY KEY,
      gene_id VARCHAR(10) NOT NULL UNIQUE           -- Entrez Gene ID
    );")
  sqliteQuickSQL(con, sql)

  gene_id <- data.frame(entrez) ## TODO: data.frame() necessary???
  sql<- paste("INSERT INTO genes(gene_id) VALUES(?);")
  dbBeginTransaction(con)
  dbGetPreparedQuery(con, sql, gene_id)
  dbCommit(con)
  message("genes table filled")
}

## The following takes a data.frame and produces a simple table from that.  It
## expects that the 1st column of that data.frame will always be entrez gene
## IDs.  All fields are assumed to be varchars of size equal to the values in
## fieldNameLens.  TODO: The cols in data have to be named and of equal
## length.  indFields is a character vector with the names of fields that we
## want indexed.  By default only _id will be indexed.
.makeEmptySimpleTable <- function(con, table, tableFieldLines){
    sql<- paste("    CREATE TABLE IF NOT EXISTS",table," (
      _id INTEGER NOT NULL,                         -- REFERENCES genes
      ",tableFieldLines,"
      FOREIGN KEY (_id) REFERENCES genes (_id)
    );") 
    sqliteQuickSQL(con, sql)
}
.makeSimpleTable <- function(data, table, con, fieldNameLens=25,
                             indFields="_id"){
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
    sql<- paste("
    INSERT INTO ",table,"
     SELECT g._id as _id, ",selFieldLines,"
     FROM genes AS g, temp AS t
     WHERE g.gene_id=t.gene_id
     ORDER BY g._id;
     ", sep="") 
    sqliteQuickSQL(con, sql)

    ## Add index to all fields in indFields (default is all)
    for(i in seq_len(length(indFields))){
    sqliteQuickSQL(con,
        paste("CREATE INDEX IF NOT EXISTS ",
              table,"_",indFields[i],"_ind ON ",table,
              " (",indFields[i],");", sep=""))      
    }
    
    ## drop the temp table
    sqliteQuickSQL(con, "DROP TABLE temp;")
  }
  message(paste(table,"table filled"))
}

## used to gather ancestor nodes for GO terms
.expandGOFrame <- function(frame, AncestMap){
  if(dim(frame)[1] ==0){
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
.filterGOFrame <- function(frame){
  message("Dropping GO IDs that are too new for the current GO.db")
  frame[frame[["go_id"]] %in% Lkeys(GOTERM),]
}

## TODO: modify this so that it no longer unwinds lists...
## used to make the 6 custom GO tables
.makeGOTablesFromNCBI <- function(con){
  bp <- .filterGOFrame(sqliteQuickSQL(con,
           paste("SELECT distinct gene_id, go_id, evidence FROM gene2go",
                 "WHERE category = 'Process'")))

  mf <- .filterGOFrame(sqliteQuickSQL(con,
           paste("SELECT distinct gene_id, go_id, evidence FROM gene2go",
                 "WHERE category = 'Function'")))
    
  cc <- .filterGOFrame(sqliteQuickSQL(con,
           paste("SELECT distinct gene_id, go_id, evidence FROM gene2go",
                 "WHERE category = 'Component'"))) 

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
    bp_all <- .expandGOFrame(bp, GOBPANCESTOR)
    mf_all <- .expandGOFrame(mf, GOMFANCESTOR)
    cc_all <- .expandGOFrame(cc, GOCCANCESTOR)
    
    .makeSimpleTable(bp_all, table = "go_bp_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))
    
    .makeSimpleTable(mf_all, table = "go_mf_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))
    
    .makeSimpleTable(cc_all, table = "go_cc_all", con, fieldNameLens=c(10,3),
                     indFields = c("_id", "go_id"))
}


#########################################################################
#########################################################################
.downloadAndSaveToTemp <- function(url, tmp){
  require(RCurl)
  if (!url.exists(url))
    stop("URL '", url, "' does not exist")
  binRes <- getBinaryURL(url)
  writeBin(binRes, con=tmp)
}

.tryDL <- function(url, tmp){
    times = 4 ## try this many times to DL
  for(i in 1:times){ 
    ## tryResult <- try( download.file(url, tmp, quiet=TRUE) , silent=TRUE)     
     tryResult <- try( .downloadAndSaveToTemp(url, tmp) , silent=TRUE)     
     if(is(tryResult,"try-error") && i < times){
         Sys.sleep(20)
     }else if(is(tryResult,"try-error") && i >= times){
         msg = paste("After 3 attempts, AnnotationDbi is still not able",
                     "to access the following URL:", url,
                     "You might want to try again later.",
                      sep=" ")
         stop(paste(strwrap(msg,exdent=2), collapse="\n"))
     }else{ return() }
  }
}



## helper for writing data to a temp NCBI database for rapid retrieval
.writeToNCBIDB <- function(NCBIcon, tableName, tmp, file){
    colClasses <- c(rep("character", times= length(unlist(file))))
    vals <- read.delim(tmp, header=FALSE, sep="\t", skip=1,
                       stringsAsFactors=FALSE, quote="",
                       colClasses = colClasses)
    if('tax_id' %in% colnames(vals)){
        vals[['tax_id']] <- as.integer(vals[['tax_id']])
    }
    colnames(vals) <- unlist(file)
    dbWriteTable(NCBIcon, name=tableName, value=vals, overwrite=TRUE)
}

.indexTaxIds <- function(NCBIcon, tableName){
    ## create INDEX on tax_id
    sql <- paste0("CREATE INDEX IF NOT EXISTS ",tableName,"_taxId ON ",
                  tableName, " (tax_id);")
    sqliteQuickSQL(NCBIcon, sql)
}


## Helper for deciding what to do about saves (returns location of saved subset)
.getFiles <- function(NCBIFilesDir, file, url, NCBIcon, tableName){
  if(is.null(NCBIFilesDir)){ ## the case for using tempfile()
      tmp <- tempfile() 
      .tryDL(url, tmp)
  }else{
      tmp <- paste(NCBIFilesDir, names(file), sep=.Platform$file.sep)
       if(!file.exists(tmp)){ 
           .tryDL(url, tmp)
       }## otherwise means we already have the file saved (AND the DB created)
  }
  .setNCBIDateStamp(NCBIcon, tableName)
  .writeToNCBIDB(NCBIcon, tableName, tmp, file)
  if("tax_id" %in% unlist(file)){
      .indexTaxIds(NCBIcon, tableName)
  }
  return(tmp)
}

## helpers for checking if NCBI 'cache database' is full (or not).
.setNCBIDateStamp <- function(NCBIcon, tableName){
    tblNm <- paste0(tableName,'_date') 
    vals = data.frame('date'=as.character(Sys.Date()))
    dbWriteTable(NCBIcon, name=tblNm, value=vals, overwrite=TRUE)
}

.getNCBIDateStamp <- function(NCBIcon, tableName){
    tblNm <- paste0(tableName,'_date')
    sqliteQuickSQL(NCBIcon, paste0("SELECT date from ",tblNm))$date
}

.isNCBICurrentWith <- function(NCBIcon, tableName){
    allTables <- dbListTables(NCBIcon)
    if(tableName %in% allTables){
        DBdate <- as.character(.getNCBIDateStamp(NCBIcon, tableName))
        curdate <- as.character(Sys.Date())
        if(DBdate == curdate){return(TRUE)}else{return(FALSE)}
    }else{return(FALSE)}
}

## helper for checking if NCBI 'cache database' is populated yet
.isNCBIPopulatedWith <- function(NCBIcon, tableName){
    exists <- dbListTables(NCBIcon)
    if(tableName %in% exists){return(TRUE)}else{return(FALSE)} 
}

##########################################################
.downloadData <- function (file, tax_id, NCBIFilesDir) {
  ## names(file) is something like: "gene2go.gz"
  message(paste("Getting data for ",names(file),sep=""))
  tableName <- sub(".gz","",names(file))
  ## Where to DL from
  url <- paste("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/",names(file), sep="")
  ## make an NCBI connection
  if(is.null(NCBIFilesDir)){
      NCBIcon <- dbConnect(SQLite(), dbname = tempfile())
  }else{
      NCBIcon <- dbConnect(SQLite(), dbname = "NCBI.sqlite")
  }
  ## if DB table is not fresh OR if table is not populated
  if(!.isNCBICurrentWith(NCBIcon, tableName) ||
     !.isNCBIPopulatedWith(NCBIcon, tableName)){
      tmp <- .getFiles(NCBIFilesDir, file, url, NCBIcon, tableName)
  }
  ## get the data from the DB
  message(paste0("extracting only data for our organism from : " ,tableName))
  if("tax_id" %in% unlist(file)){ ## when there is a tax_id need to subset
      sql <- paste0("SELECT * FROM ", tableName, " WHERE ",
                    tableName, ".tax_id = '", tax_id, "'")
      vals <- sqliteQuickSQL(NCBIcon, sql)
      ## if things are empty: make sure we have correct number of cols
      if(dim(vals)[1] == 0){## if after select there are no records here
          vals <- data.frame(t(1:length(unlist(file))),
                             stringsAsFactors=FALSE)
          vals <- vals[FALSE,]
          colnames(vals) <- unlist(file)
      }
  }else{ ## Just get the whole record
      message(paste0("getting all data for our organism from : " ,tableName))
      sql <- paste0("SELECT * FROM ", tableName)
      vals <- sqliteQuickSQL(NCBIcon, sql)
  }
  ## remove row_names col 
  vals[,!(colnames(vals) %in% 'row_names')]
}




## helper to populate base tables
.populateBaseTable <- function(con, sql, data, table) {
    dbBeginTransaction(con)
    dbGetPreparedQuery(con, sql, data)
    dbCommit(con)
    message(paste("table ",table," filled",sep=""))
}

## Convenience function for quickly making indices
.makeIndex <- function(con,table, field){
  indName <-paste("ind",table,field,sep="_")
  sqliteQuickSQL(con,
   paste("CREATE INDEX IF NOT EXISTS",indName,"ON",table,"(",field,")"))
}

## metadata stuff
.createMetadataTables <- function(con){
  sqliteQuickSQL(con, paste("CREATE TABLE IF NOT EXISTS metadata (name",
                            "VARCHAR(80) PRIMARY KEY,value VARCHAR(255))"))  
  sqliteQuickSQL(con, paste("CREATE TABLE IF NOT EXISTS map_metadata ",
                            "(map_name VARCHAR(80) NOT NULL,",
                            "source_name VARCHAR(80) NOT NULL,",
                            "source_url VARCHAR(255) NOT NULL,",
                            "source_date VARCHAR(20) NOT NULL)"))
  sqliteQuickSQL(con, paste("CREATE TABLE IF NOT EXISTS map_counts ",
                            "(map_name VARCHAR(80) PRIMARY KEY,",
                            "count INTEGER NOT NULL)"))    
}

.addMeta <- function(con, name, value){
  data = data.frame(name,value,stringsAsFactors=FALSE)
  sql <- "INSERT INTO metadata (name,value) VALUES(?,?)"
  .populateBaseTable(con, sql, data, "metadata")
}

.addMetadata <- function(con, tax_id, genus, species){
  name <- c("DBSCHEMAVERSION","DBSCHEMA","ORGANISM","SPECIES","CENTRALID",
            "TAXID",
            "EGSOURCEDATE","EGSOURCENAME","EGSOURCEURL",
            "GOSOURCEDATE","GOSOURCENAME","GOSOURCEURL",
            "GOEGSOURCEDATE","GOEGSOURCENAME","GOEGSOURCEURL",
            "Db type","Supporting package")
  value<- c("2.1","ORGANISM_DB",paste(genus,species),paste(genus,species),"EG",
            tax_id,
       date(), "Entrez Gene","ftp://ftp.ncbi.nlm.nih.gov/gene/DATA",
       .getGODate(), "Gene Ontology","ftp://ftp.geneontology.org/pub/go/godata",
       date(), "Entrez Gene","ftp://ftp.ncbi.nlm.nih.gov/gene/DATA",
            "OrgDb","AnnotationDbi")
  .addMeta(con, name, value)
}

.addMapMeta <- function(con, map_name, source_name, source_url, source_date){
  data = data.frame(map_name,source_name,source_url,
    source_date,stringsAsFactors=FALSE)
  sql <- paste("INSERT INTO map_metadata (map_name,source_name,source_url,",
               "source_date) VALUES(?,?,?,?)")
  .populateBaseTable(con, sql, data, "map_metadata")
}

.addMapMetadata <- function(con, tax_id, genus, species){
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
  goSrcs <- sqliteQuickSQL(con,
             "SELECT source_name FROM map_metadata WHERE map_name='GO2EG'")
  if("Gene Ontology" %in% goSrcs & "blast2GO" %in% goSrcs){
    sqliteQuickSQL(con,paste("DELETE FROM map_metadata where map_name",
                             "='GO2EG' AND source_name='Gene Ontology'"))
  }
}


## Use when there is no GO data
.getBlast2GOData <- function(tax_id, con) {
  url = paste("http://www.b2gfar.org/_media/species:data:",
        tax_id,".annot.zip",sep="")
  tmp <- tempfile()
  ##download.file(url, tmp, quiet=TRUE)
  .tryDL(url,tmp)
  vals <- read.delim(unzip(tmp), header=FALSE, sep="\t", quote="",
                     stringsAsFactors=FALSE)
  ## I will need to so extra stuff here to match up categories etc.
  ## (vals has to look like gene2go would, and I have to join to refseq and to
  ## accession just to get my EGs)
  RSVals <- vals[grep("refseq",vals[,4]),c(2,3)]
  GBVals <- vals[grep("genbank",vals[,4]),c(2,6)]
  colnames(GBVals) <- colnames(RSVals) <- c("go_id","accession")
  tables <- sqliteQuickSQL(con,"SELECT * FROM sqlite_master")$name

  if(!any(c("gene2refseq","gene2accession") %in% tables)){
    stop("It is impossible to match up blasted GO terms with NCBI accessions.")
  }else{
      if("gene2refseq" %in% tables){
          g2rs <- sqliteQuickSQL(con,
                                 "SELECT gene_id, protein_accession FROM gene2refseq")
          g2rs[,2] <- sub("\\.\\d+$","",g2rs[,2])
          vals <- merge(RSVals,g2rs, by.x="accession" , by.y="protein_accession")
      }
      if("gene2accession" %in% tables){
          g2ac <- sqliteQuickSQL(con,
                                 "SELECT gene_id, protein_accession FROM gene2accession")
          gb <- merge(GBVals,g2ac, by.x="accession" , by.y="protein_accession")
          if("gene2refseq" %in% tables){
              vals <- rbind(vals, gb)
          }else{
              vals <- gb
          }
      }
  }

  ## Add to the metadata:
  name <- c("BL2GOSOURCEDATE","BL2GOSOURCENAME","BL2GOSOURCEURL")
  value<- c(date(),"blast2GO","http://www.blast2go.de/")
  .addMeta(con, name, value)

  ## modify and then add back to map_metadata
  map_name <- c("GO2EG", "GO2ALLEGS")
  source_name <- c(rep("blast2GO",2))
  source_url <- c(rep("http://www.blast2go.de/",2))
  source_date <- c(rep(date(),2))
  .addMapMeta(con, map_name, source_name, source_url, source_date)
    
  ## then assemble the data back together to make it look like gene2go so
  ## it can be inserted there.
  rlen <- dim(vals)[1]
  data <-data.frame(tax_id = rep(tax_id, rlen),
                   gene_id = as.character(vals[["gene_id"]]),
                   go_id = vals[["go_id"]],
                   evidence = rep("IEA", rlen),
                   go_qualifier = rep("-", rlen),
                   go_description = rep("-", rlen),
                   pubmed_id = rep("-", rlen),
                   category = Ontology(vals[["go_id"]]),
                   stringsAsFactors=FALSE)
  data[,8] <- sub("BP","Process",data[,8])
  data[,8] <- sub("CC","Component",data[,8])
  data[,8] <- sub("MF","Function",data[,8])
  ## cleanup of file left behind by unzip() operation.
  if(file.exists(paste(tax_id,".annot",sep=""))){
    file.remove(paste(tax_id,".annot",sep=""))
  }
  ## return data
  data
}

## need to be able to generate the columns data from the files and the types.
.generateCols <- function (file) {
  types <- c("tax_id" = "INTEGER NOT NULL",
             "gene_id" = "INTEGER NOT NULL",
             "go_id" = "TEXT NOT NULL",
             "evidence" = "TEXT",
             "go_qualifier" = "TEXT",
             "go_description" = "TEXT",
             "pubmed_id" = "INTEGER",
             "category" = "TEXT",
             "status" = "TEXT",
             "rna_accession" = "TEXT",
             "rna_gi" = "INTEGER",
             "protein_accession" = "TEXT",
             "protein_gi" = "INTEGER",
             "genomic_dna_accession" = "TEXT",
             "genomic_dna_gi" = "INTEGER",
             "genomic_start" = "INTEGER",
             "genomic_end" = "INTEGER",
             "orientation" = "TEXT",
             "assembly" = "TEXT",
             "unigene_id" = "TEXT",
             "symbol" = "TEXT",
             "locus_tag" = "TEXT",
             "synonyms" = "TEXT",
             "dbXrefs" = "TEXT",
             "chromosome" = "TEXT",
             "map_location" = "TEXT",
             "description" = "TEXT",
             "gene_type" = "TEXT",
             "nomenclature_symbol" = "TEXT",
             "nomenclature_name" = "TEXT",
             "nomenclature_status" = "TEXT",
             "other_designations" = "TEXT",
             "mim_id" = "INTEGER NOT NULL",
             "relation_type" = "TEXT",
             "refseq_id" = "TEXT NOT NULL",
             "uniprot_id" = "TEXT NOT NULL"
             )
  cols <- character()
  for(i in seq_len(length(file[[1]]))){
    cols[i] <-  paste(file[[1]][i],types[file[[1]]][i])
  }
  paste(cols, collapse =",")
}

.mergeAndCleanAccession <- function(ac){
  ## pair 1st col with 2nd, clean and rename
  acr <- ac[ac[,2]!="-",1:2]
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
.generateTempTableINSERTStatement <- function (file) {
  Qs <- paste(rep("?",length(file[[1]])),  collapse=",")
  paste("INSERT INTO ",sub(".gz$","",names(file)),
  " VALUES(",Qs,");",sep="")
}


## need code to generate a table for each of the
## file below is like: files[1] or files[2]
.createTEMPNCBIBaseTable <- function (con, file, tax_id, NCBIFilesDir) {
  data <- .downloadData(file, tax_id, NCBIFilesDir=NCBIFilesDir)
  table <- sub(".gz$","",names(file))
  if(is.null(table)) stop("Unable to infer a table name.")
  cols <- .generateCols(file)
  sql<- paste("CREATE TABLE IF NOT EXISTS ",table," (",cols,");")
  sqliteQuickSQL(con, sql)
  
  message(paste("Populating ",table," table:", sep=""))
  sql <- .generateTempTableINSERTStatement(file)
  if(dim(data)[1] != 0){ ## ie. there has to be SOME data...
    .populateBaseTable(con, sql, data, table)
  }else if(dim(data)[1] == 0 && names(file)=="gene2go.gz"){ ## we need blast2GO
    message(paste("Getting blast2GO data as a substitute for ",table,sep=""))
    data <- .getBlast2GOData(tax_id, con)
    sql <- paste("INSERT INTO gene2go(tax_id, gene_id, go_id, evidence, ",
                 "go_qualifier, go_description,",
                 "pubmed_id, category) VALUES(?,?,?,?,?,?,?,?);")
    .populateBaseTable(con, sql, data, table)
  }else{
    message(paste("No data available for table ",table,sep=""))
  }
}

## Download all the data files and make the cache DB
.setupBaseDBFromDLs <- function(files, tax_id, con, NCBIFilesDir){
  for(i in seq_len(length(files))){
    .createTEMPNCBIBaseTable(con, files[i], tax_id, NCBIFilesDir=NCBIFilesDir)
  }
  con
}



.dropOldTables <- function(con,fileNames){
  fileNames <- sub(".gz", "", fileNames)
  message(paste("dropping table",fileNames))
  for(i in seq_len(length(fileNames))){
    sqliteQuickSQL(con, paste("DROP TABLE IF EXISTS",fileNames[i]))
  }
}


.generateOrgDbName <- function(genus, species){
    specStr <- paste(toupper(substr(genus,1,1)),species,sep="")
    paste("org.",specStr,".eg",sep="")
}

.getGODate <- function(){
  GOinf <- dbInfo(GO_dbconn())
  GOinf[GOinf[["name"]]=="GOSOURCEDATE","value"]
}


## wanted when we want to know how many of something we can get
.computeSimpleEGMapCounts <- function(con, table, field){
  sql<- paste("SELECT count(DISTINCT g.gene_id) FROM ",table,
              " AS t, genes as g WHERE t._id=g._id AND t.",field,
              " NOT NULL", sep="")
  message(sql)
  sqliteQuickSQL(con,sql)
}

## when we want to know how many of something there is (reverse mapping?)
.computeSimpleMapCounts <- function(con, table, field){
  sql<- paste("SELECT count(DISTINCT t.",field,") FROM ",table,
              " AS t, genes as g WHERE t._id=g._id AND t.",field,
              " NOT NULL", sep="")
  message(sql)
  sqliteQuickSQL(con,sql)
}


## TODO: correct the counts as needed
.addMapCounts <- function(con, tax_id, genus, species){
  map_name <- c("GENENAME","SYMBOL","SYMBOL2EG","CHR","REFSEQ","REFSEQ2EG",
                "UNIGENE","UNIGENE2EG","GO","GO2EG","GO2ALLEGS","ALIAS2EG",
                "TOTAL")
  count <- c(.computeSimpleEGMapCounts(con, "gene_info", "gene_name"),
             .computeSimpleEGMapCounts(con, "gene_info", "symbol"),
             .computeSimpleMapCounts(con, "gene_info", "symbol"),
             .computeSimpleEGMapCounts(con, "chromosomes", "chromosome"),
             .computeSimpleEGMapCounts(con, "refseq", "accession"),
             .computeSimpleMapCounts(con, "refseq", "accession"),
             .computeSimpleEGMapCounts(con, "unigene", "unigene_id"),
             .computeSimpleMapCounts(con, "unigene", "unigene_id"),
             .computeSimpleEGMapCounts(con, "accessions", "accession"),
             .computeSimpleMapCounts(con, "accessions", "accession"),
             #GO ones 
             .computeSimpleEGMapCounts(con, "alias", "alias_symbol"),
             sqliteQuickSQL(con,"SELECT count(DISTINCT gene_id) FROM genes")
             )
  data = data.frame(map_name,count,stringsAsFactors=FALSE)
  sql <- "INSERT INTO map_counts (map_name,count) VALUES(?,?)"
  .populateBaseTable(con, sql, data, "map_counts")
}

## list of primary NCBI files that we need to process
.primaryFiles <- function(){
    list(
         "gene2pubmed.gz" = c("tax_id","gene_id", "pubmed_id"),
         "gene2accession.gz" = c("tax_id","gene_id","status","rna_accession",
           "rna_gi","protein_accession","protein_gi","genomic_dna_accession",
           "genomic_dna_gi","genomic_start","genomic_end","orientation",
           "assembly"),
         ## This one might be needed later
         "gene2refseq.gz" = c("tax_id","gene_id","status","rna_accession",
           "rna_gi","protein_accession","protein_gi","genomic_dna_accession",
           "genomic_dna_gi","genomic_start","genomic_end","orientation",
           "assembly"),
         "gene2unigene" = c("gene_id","unigene_id"),
         "gene_info.gz" = c("tax_id","gene_id","symbol","locus_tag",
           "synonyms","dbXrefs","chromosome","map_location","description",
           "gene_type","nomenclature_symbol","nomenclature_name",
           "nomenclature_status","other_designations", "modification_date"),
         ##        "mim2gene.gz" = c("mim_id","gene_id","relation_type"),
         ##        "gene_refseq_uniprotkb_collab.gz" =
         ##         c("refseq_id","uniprot_id"),
         "gene2go.gz" = c("tax_id","gene_id","go_id","evidence",
           "go_qualifier", "go_description","pubmed_id","category")
         )
}


#########################################################################
## Generate the database using the helper functions:
#########################################################################

makeOrgDbFromNCBI <- function(tax_id, genus, species, NCBIFilesDir,
                              outputDir){
  require(RSQLite)
  require(GO.db)
  dbFileName <- paste0(.generateOrgDbName(genus,species),".sqlite")
  dbFileName <- file.path(outputDir, dbFileName)
  if(file.exists(dbFileName)){ file.remove(dbFileName) }
  con <- dbConnect(SQLite(), dbFileName)
  .createMetadataTables(con)  ## just makes the tables
  ## I need a list of files, along with their column names 
  ## (needed for schema definitions later)
  ## IF ANY OF THESE gets moved in the source files
  ## then things will be in the wrong place!
  files = .primaryFiles()
  .setupBaseDBFromDLs(files, tax_id, con, NCBIFilesDir=NCBIFilesDir)
  
  ## Add metadata:
  .addMetadata(con, tax_id, genus, species) 
  ## Add map_metadata:
  .addMapMetadata(con, tax_id, genus, species)
 
  ## Make the central table
  egs <- sqliteQuickSQL(con, "SELECT distinct gene_id FROM gene_info")[,1]
  .makeCentralTable(egs, con)
  
  ## Make the other tables:
  ## gene_info
  symbs <- sqliteQuickSQL(con,
    "SELECT distinct gene_id, description, symbol FROM gene_info")
  colnames(symbs) <- c("gene_id", "gene_name", "symbol")
  ## constraint requires that we always have BOTH a symbol and a name.
  symbs <- symbs[!is.na(symbs[,2]) & !is.na(symbs[,3]),] # strict constraint!
  ## Only 10 things fail this strict constraint for human.
  ## ALL because of no gene symbol (but that have a name)
  ## TODO: fix this so that I can make this table without
  ## such a strict constraint
  .makeSimpleTable(symbs, table="gene_info_temp", con, fieldNameLens=c(255,80))
  
  ## alias ## requires sub-parsing.
  alias <- sqliteQuickSQL(con,
    "SELECT distinct gene_id, synonyms FROM gene_info")
  aliases <- sapply(alias[,2], strsplit, "\\|")
  numAlias <- sapply(aliases, length)
  alGenes <- rep(alias[,1], numAlias)
  alias <- data.frame(gene_id=alGenes,alias_symbol=unlist(aliases))
  .makeSimpleTable(alias, table="alias", con)
  
  ## chr
  chrs <- sqliteQuickSQL(con,
    "SELECT distinct gene_id, chromosome FROM gene_info")
    .makeSimpleTable(chrs, table="chromosomes", con)
  
  
  ## pubmed
  pm <- sqliteQuickSQL(con,
    "SELECT distinct gene_id, pubmed_id FROM gene2pubmed")
    .makeSimpleTable(pm, table="pubmed", con)
  
  ## refseq ## requires sub-parsing.
  rs <- sqliteQuickSQL(con,
    "SELECT distinct gene_id,rna_accession,protein_accession FROM gene2refseq")
  rs <- .mergeAndCleanAccession(rs)
  .makeSimpleTable(rs, table="refseq", con)
  
  ## accessions
  ac <- sqliteQuickSQL(con,
     paste("SELECT distinct gene_id,rna_accession,protein_accession",
           "FROM gene2accession"))
  ac <- .mergeAndCleanAccession(ac)
  .makeSimpleTable(ac, table="accessions", con)
  
  ## unigene 
  ug <- sqliteQuickSQL(con,
    "SELECT distinct g.gene_id, u.unigene_id FROM gene2unigene as u,
     gene_info as g WHERE u.gene_id=g.gene_id")
  .makeSimpleTable(ug, table="unigene", con) 
  
  ## Make the GO tables:
  .makeGOTablesFromNCBI(con) 
  
  ## Drop all the older tables (which will include the original "gene_info").
  .dropOldTables(con,names(files))  
  ## Rename "gene_info_temp" to be just "gene_info":
  sqliteQuickSQL(con,"ALTER TABLE gene_info_temp RENAME TO gene_info")
  ## add GO views to the DB
  makeGOViews(con)

  ## Add map_counts: (must happen at end)
  .addMapCounts(con, tax_id, genus, species)
  
}



## For argument checking:
.isSingleString <- function(x){
    is.atomic(x) && length(x) == 1L && is.character(x)
}
.isSingleStringOrNA <- function(x)
{
    is.atomic(x) && length(x) == 1L && (is.character(x) || is.na(x))
}
.isSingleStringOrNull <- function(x)
{
    (is.atomic(x) && length(x) == 1L && is.character(x)) ||
    (is.atomic(x) && length(x) == 0L && is.null(x))
}



## OLDER function to make the package:
OLD_makeOrgPackageFromNCBI <- function(version,
                               maintainer,
                               author,
                               outputDir,
                               tax_id,
                               genus,
                               species,
                               NCBIFilesDir){
  message("If this is the 1st time you have run this function, it may take a long time (over an hour) to download needed files and assemble a 12 GB cache databse in the NCBIFilesDir directory.  Subsequent calls to this function should be faster (seconds).  The cache will try to rebuild once per day.")
  ## Arguement checking:
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
  if(!.isSingleStringOrNull(NCBIFilesDir))
      stop("'NCBIFilesDir' argument needs to be a single string or NULL")
       
  makeOrgDbFromNCBI(tax_id=tax_id, genus=genus, species=species,
                    NCBIFilesDir=NCBIFilesDir, outputDir)
  
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



#############################################################################
## Code specific for the new process.

## loop through and make the tables
.makeBaseDBFromDLs <- function(files, tax_id, con, NCBIFilesDir){
    data <- list()
    for(i in seq_len(length(files))){
        res <- .downloadData(files[i], tax_id, NCBIFilesDir=NCBIFilesDir)
        data[[i]] <- res
    }
    data
}


## Helper functions for extracting data from a frame, making sure
## there *is* some data and renaming it etc.
.tossHyphens <- function(df){
    cols <- 1:dim(df)[2] ## get range of columns in df
    cols <- cols[-1] ## drop the 1st one
    for(i in cols){
        df <- df[df[[i]]!='-',]
    }
    df
}

## extraction function
.extract <- function(rawData, keepTable, keepCols){
    rawData[[keepTable]][,keepCols]
}

## rename and alter the data list (if needed)
.rename <- function(data, dataSub, newName, newCols){
    dataSub <- .tossHyphens(dataSub)
    if(dim(dataSub)[1] > 0){ ## if there are rows...
        data[[newName]] <- unique(dataSub)
        colnames(data[[newName]]) <- newCols
        return(data)
    }else{ ## don't add anything to data
        return(data)
    }
}

## extract AND rename etc.
.procDat <- function(rawData, keepTable, keepCols, data, newName, newCols){
    pubSub <- .extract(rawData, keepTable, keepCols)
    .rename(data, pubSub, newName, newCols)
}


## based (very loosely) on older .Blast2GOData() above
## Alternative Blast2GO (refactor in progress)
## The file we have to proces has specs here:
## http://www.b2gfar.org/fileformat

.getBlast2GO <- function(tax_id, refseq, accs) {
  url = paste("http://www.b2gfar.org/_media/species:data:",
        tax_id,".annot.zip",sep="")
  tmp <- tempfile()
  ##download.file(url, tmp, quiet=TRUE)
  .tryDL(url,tmp)
  rawVals <- read.delim(unzip(tmp), header=FALSE, sep="\t", quote="",
                     stringsAsFactors=FALSE)
  ## I will need to so extra stuff here to match up categories etc.
  ## (vals has to look like gene2go would, and I have to join to refseq and to
  ## accession just to get my EGs)
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


## Helper to get complete list of entrez gene IDs
.mkGIdQ <- function(tax_id, tableName){
    paste0("SELECT gene_id FROM ", tableName, " WHERE ",
           tableName, ".tax_id = '", tax_id, "'")
}
.asV <- function(df){
    unique(as.character(t(df)))
}
.getAllEnrezGeneIdsFromNCBI <- function(tax_id, NCBIcon){
    id1 <- .asV(sqliteQuickSQL(NCBIcon, .mkGIdQ(tax_id, 'gene_info')))
    id2 <- .asV(sqliteQuickSQL(NCBIcon, .mkGIdQ(tax_id, 'gene2refseq')))
    id3 <- .asV(sqliteQuickSQL(NCBIcon, .mkGIdQ(tax_id, 'gene2go')))
    id4 <- .asV(sqliteQuickSQL(NCBIcon, .mkGIdQ(tax_id, 'gene2pubmed')))
    id5 <- .asV(sqliteQuickSQL(NCBIcon, .mkGIdQ(tax_id, 'gene2accession')))
    unique(c(id1,id2,id3,id4,id5))
}


## Helper to make a list of data frames from the NCBI cache database
## (which it will also make if needed)
prepareDataFromNCBI <- function(tax_id=tax_id, NCBIFilesDir=NCBIFilesDir,
                                outputDir){
    ## Get list of files that need to be processed 
    files = .primaryFiles()
    NCBIcon <- dbConnect(SQLite(), dbname = "NCBI.sqlite")
    rawData <- .makeBaseDBFromDLs(files, tax_id, NCBIcon, NCBIFilesDir)
    names(rawData) <- names(files)
    data = list() ## so I can pass them back
    ## now post-process the data list into proper data frames for each
    
    ## pubmed
    data <- .procDat(rawData,
                     keepTable='gene2pubmed.gz',
                     keepCols=c('gene_id','pubmed_id'),
                     data,
                     newName='pubmed',
                     newCols=c('GID', 'PMID'))
    ## chromosomes
    data <- .procDat(rawData,
                     keepTable='gene_info.gz',
                     keepCols=c('gene_id','chromosome'),
                     data,
                     newName='chromosomes',
                     newCols=c('GID', 'CHR'))
    ## gene_info
    data <- .procDat(rawData,
                     keepTable='gene_info.gz',
                     keepCols=c('gene_id','description','symbol'),
                     data,
                     newName='gene_info',
                     newCols=c('GID', 'GENENAME','SYMBOL'))
    
    ## entrez genes (special code to ensure we have them all)
    egIDs <- .getAllEnrezGeneIdsFromNCBI(tax_id, NCBIcon)
    egData <- data.frame(GID=egIDs, ENTREZID=egIDs,stringsAsFactors=FALSE)
    data <- .rename(data,
                    egData,
                    newName='entrez_genes',
                    newCols=c('GID', 'ENTREZID'))
    
    ## Alias table needs to contain both symbols AND stuff from synonyms
    ## For alias I need to massage the synonyms field from gene_info
    ## and add it to symbols from the same.
    aliasSQL <- paste0("SELECT distinct gene_id, synonyms FROM gene_info ",
                       "WHERE tax_id ='",tax_id,"'")
    alias <- sqliteQuickSQL(NCBIcon,aliasSQL)
    aliases <- sapply(alias[,2], strsplit, "\\|")
    numAlias <- sapply(aliases, length)
    allGenes <- rep(alias[,1], numAlias)
    aliasDat <- data.frame(gene_id=allGenes,alias_symbol=unlist(aliases)
                           ,stringsAsFactors=FALSE)
    symbolDat <- .extract(rawData,
                          keepTable='gene_info.gz',
                          keepCols=c('gene_id','symbol'))
    colnames(aliasDat) <- NULL
    colnames(symbolDat) <- NULL
    faliasDat <- unique(rbind(as.matrix(aliasDat), as.matrix(symbolDat)))
    rownames(faliasDat) <- NULL
    data <- .rename(data,
                    as.data.frame(faliasDat, stringsAsFactors=FALSE),
                    newName='alias',
                    newCols=c('GID', 'ALIAS'))            
    ## refseq requires a custom job since a couple columns must be
    ## combined into one column
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
    ## accession is similar to refseq
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
    


    ## go has to be done in two separate steps so that if we get
    ## nothing back we can look for results from blast2GO.
    goDat <- .extract(rawData,
                      keepTable='gene2go.gz',
                      keepCols=c('gene_id','go_id','evidence'))
    if(dim(goDat)[1] == 0){## then try Blast2GO
        ## get all accessions together
        goDat <- .getBlast2GO(tax_id, data[['refseq']], data[['accessions']]) 
    }
    ## .rename already checks to make sure there actually ARE GO IDs
    data <- .rename(data,
                    goDat,
                    newName='go',
                    newCols=c("GID","GO","EVIDENCE"))
    
    
    ## unigene needs custom extraction of relevant genes ONLY (do last)
    uniDat <- .extract(rawData,
                       keepTable='gene2unigene',
                       keepCols=c('gene_id','unigene_id'))    
    ## Then row filter our things that are not in our set of entrez gene IDs
    uniDat <- uniDat[uniDat$gene_id %in% egIDs,]
    ## and rename
    data <- .rename(data,
                    dataSub=uniDat,
                    newName='unigene',
                    newCols=c('GID', 'UNIGENE'))
    
    
    ## TODO: Then make sure we have key things (gene names, GO IDs etc.)
    
    

    data
}


## NEW function to make the package:
NEW_makeOrgPackageFromNCBI <- function(version,
                               maintainer,
                               author,
                               outputDir,
                               tax_id,
                               genus,
                               species,
                               NCBIFilesDir){
  message("If this is the 1st time you have run this function, it may take a long time (over an hour) to download needed files and assemble a 12 GB cache databse in the NCBIFilesDir directory.  Subsequent calls to this function should be faster (seconds).  The cache will try to rebuild once per day.")
  ## Arguement checking:
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
  if(!.isSingleStringOrNull(NCBIFilesDir))
      stop("'NCBIFilesDir' argument needs to be a single string or NULL")

  data <- prepareDataFromNCBI(tax_id=tax_id, NCBIFilesDir=NCBIFilesDir,
                              outputDir)
  
  dbName <- .generateOrgDbName(genus,species)
  dbfile <- paste0(dbName, ".sqlite")

  ## if there is go data, then we have to use the goTable argument.
  if(!is.null(data[['go']])){
      .makeOrgPackage(data,
                      version=version,
                      maintainer=maintainer,
                      author=author,
                      outputDir=outputDir,
                      tax_id=tax_id,
                      genus=genus,
                      species=species,
                      goTable="go") 
  }else{
      .makeOrgPackage(data,
                      version=version,
                      maintainer=maintainer,
                      author=author,
                      outputDir=outputDir,
                      tax_id=tax_id,
                      genus=genus,
                      species=species) 
  }
}





## function wrapper to make the package:
makeOrgPackageFromNCBI <- function(version,
                                   maintainer,
                                   author,
                                   outputDir=getwd(),
                                   tax_id,
                                   genus,
                                   species,
                                   NCBIFilesDir=getwd(),
                                   useDeprecatedStyle=FALSE){
    if(useDeprecatedStyle==TRUE){
        OLD_makeOrgPackageFromNCBI(version,maintainer,author,outputDir,
                                   tax_id,genus,species,NCBIFilesDir)
    }else{
        NEW_makeOrgPackageFromNCBI(version,maintainer,author,outputDir,
                                   tax_id,genus,species,NCBIFilesDir)
    }
}




## STILL TODO:

## 4- finish code to check that we have all the parts and make sure it
## filters out any stuff that we don't need (empty data.frames should
## be dropped) = For this we want to add an argument that normally
## will be FALSE that will be something like requireMinStandards and
## which will not actually make a DB (will error out) in the case
## where certain standards are not met...

## 5- extracting ALL unigenes into memory every time is not efficient.  It would be better if I knew beforehand which organisms I need this for (and could therefore not do it the rest of the time).
