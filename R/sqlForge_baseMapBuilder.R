## Just a utility to prevent empty IDs from ever causing any more mayhem
cleanSrcMap <- function(file) {
    fileVals <- read.delim(file=file, header=FALSE, sep="\t", quote="", colClasses = "character")
    fileVals <- fileVals[fileVals[,1]!='',]

    ##For the case where someone puts no value in, we need things to be an NA
    fileVals[!is.na(fileVals[,2]) & fileVals[,2]=='',2] <- NA
    
    ##Expand IDs to match all the ones available
    probe <- fileVals[,1]
    id <- fileVals[,2]
    id <- gsub(" ","", id)
    id <- strsplit(id, split=";", fixed=TRUE)
    ##Otherwise, the following may remove probes that map to nothing
    id_count <- sapply(id, length)  
    id_probe <- rep(probe, id_count)
    id <- unlist(id)

    insVals <- cbind(id_probe, id)
    insVals <- as.data.frame(insVals)
    insVals
}

cleanRefSeqs <- function(baseMap){
    baseMap = as.matrix(baseMap)
    baseMap[,2] = sub("\\.\\d+?$", "", baseMap[,2], perl=TRUE)
    baseMap = cbind(baseMap[,1],baseMap[,2])
    baseMap = as.data.frame(baseMap)
    baseMap
}

printOutBaseMaps <- function(baseMap, pkgName, otherSrcObjs){
    write.table(baseMap, paste(pkgName,"_baseMap.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    for(i in 1:length(otherSrcObjs)){
        write.table(otherSrcObjs[i], paste(pkgName,"_otherSrcMap_",i,".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }
}


## scrubs Affy file of comments and returns a data.frame that is free of them
scrubAffyFileComments <- function(file){
  ## 1st use readLines to pull in file and clean out the cruft.
  con <- file(file)
  rows <- readLines(con)
  filtInd <- grep("^#", rows)
  if(length(filtInd)>0){## check if there were actually any comments
    rows <- rows[-filtInd] ##if so, remove those lines.
  }
  tmpFl <- tempfile()
  write(rows, file=tmpFl)
  res <- read.csv(tmpFl, as.is=TRUE, na.strings="---", colClasses="character")
  unlink(tmpFl)
  close(con)
  res
}

makeBaseMaps <- function(csvFileName,
                          GenBankIDName="Representative.Public.ID",
                          EntrezGeneIDName="Entrez.Gene",
                          targetDir="."
                          ) {
    ##csvFile <- read.csv(csvFileName, as.is=TRUE, na.strings="---", colClasses="character")
    csvFile <- scrubAffyFileComments(csvFileName)
      
    probe <- csvFile[,1]
    gb <- csvFile[, GenBankIDName]
    eg <- csvFile[, EntrezGeneIDName]
    rm(csvFile)
    gb <- lapply(unlist(gb), function(x) toupper(strsplit(x,"\\.")[[1]][1]))
    gb <- unlist(gb)
    eg <- strsplit(eg, split=" /// ", fixed=TRUE)
    eg_count <- sapply(eg, length)
    eg_probe <- rep(probe, eg_count)
    eg <- unlist(eg)
    
    baseFiles = list()
    
    baseFiles[[1]] = cbind(probe, gb)
    baseFiles[[2]] = cbind(eg_probe, eg)
    baseFiles
}


##For labelDuplicatedProbes() The 1st col HAS to be probes, and the 2nd col
##MUST be the central ID for the platform!
labelDuplicatedProbes <- function(frame){
    ##I need to test 1st if the probe matches more than one thing in the 1st
    ##col (numProbes) Then I need to test if among those matches there is more
    ##than one entrez gene ID.  If more than one, then I need to put a 1
    ##There is a corner case however, if I have 1 unmapped and one
   tmp <- tapply(1:nrow(frame), frame[,1], function(x) frame[x,2])
   tmp2 <- sapply(tmp, length)
   tmp <- sapply(tmp, function(x) length(unique(x[!is.na(x)])))
   out <- tmp > 1 & tmp2 > 1
   out <- out[frame[,1]]
   is_multiple <- as.numeric(out)
   out <- cbind(frame, is_multiple)
   out
}


probe2gene <- function(baseMap, otherSrc,
			baseMapType=c("gb", "ug", "eg", "refseq", "gbNRef", "image"), 
			chipMapSrc, chipSrc, pkgName, outputDir=".", allowHomologyGene=FALSE) {
        ## message(cat(paste("Using '",chipSrc,"' for chipSrc.", sep="")))
	baseMapType <- match.arg(baseMapType)
        require("RSQLite")
	drv <- dbDriver("SQLite")
	outputFile <- file.path(outputDir, paste(pkgName, "sqlite", sep="."))	
	db <- dbConnect(drv, outputFile)
	dbGetQuery(db, 
		"CREATE TABLE metadata (name VARCHAR(80) PRIMARY KEY, value VARCHAR(255) );")
	metadata_sql <- paste("INSERT INTO metadata VALUES ('PKGNAME', '",
				pkgName, "');", sep="", collapse="")
	dbGetQuery(db, metadata_sql)
	dbGetQuery(db, "CREATE TABLE curr_map (probe_id TEXT, gene_id TEXT);")  #curr_map is the baseMap file
	dbGetQuery(db, "CREATE TABLE probes_ori (probe_id TEXT);")             #all the probes from curr_map
	dbGetQuery(db, "CREATE TABLE other (probe_id TEXT, gene_id TEXT);")     #This holds the otherSrc data
	dbGetQuery(db, "CREATE TABLE probe2gene (probe_id TEXT, gene_id TEXT);") 
	dbGetQuery(db, "CREATE TABLE probe2acc (probe_id TEXT, gene_id TEXT);") 
	dbGetQuery(db, "CREATE TABLE temp_probe_map (probe_id TEXT, gene_id TEXT, accession TEXT);")
	dbGetQuery(db, "CREATE TABLE probe_map (probe_id TEXT, gene_id TEXT, accession TEXT, is_multiple SMALLINT NOT NULL);")

        #If there might be refseq IDs, then they might have unwanted .<digit> extensions which must be removed.
        #This (unfortunately cannot be done for the otherSrc files since I won't know what I am getting in that case).
        if(baseMapType=='refseq' || baseMapType=='gbNRef'){  #need to verify that this function is not destructive to genbank IDs
            baseMap=cleanRefSeqs(baseMap)
        }
        
        #populate the contents of baseMap into curr_map
        clnVals <-baseMap
        sqlIns <- "INSERT INTO curr_map (probe_id,gene_id) VALUES (?,?);"
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sqlIns, clnVals)
        dbClearResult(rset)
        dbCommit(db)

        dbGetQuery(db, "INSERT INTO probes_ori SELECT DISTINCT probe_id FROM curr_map;")
	dbGetQuery(db,
		"DELETE FROM curr_map WHERE gene_id='NA' OR gene_id='';")
	attach_sql <- paste("ATTACH DATABASE '",chipMapSrc,"' AS src;", sep="", collapse="")
	dbGetQuery(db, attach_sql)
	if(baseMapType=='eg') {
                message(cat("baseMapType is eg"))
		dbGetQuery(db, "INSERT INTO probe2acc SELECT probe_id, NULL FROM probes_ori;")
		sql <- paste("INSERT INTO probe2gene",
                             "SELECT probe_id, gene_id",  ##To make things map to multiples, I will have to remove the DISTINCT and GROUP BY clauses when I insert into probe2gene
                             "FROM curr_map;", sep=" ", collapse="")
		dbGetQuery(db, sql)
	} else if (baseMapType=='image') {
                message(cat("baseMapType is image"))
		dbGetQuery(db, "INSERT INTO probe2acc SELECT probe_id, NULL FROM probes_ori;")
		sql <- paste("INSERT INTO probe2gene",
			    "SELECT c.probe_id, i.gene_id",
			    "FROM curr_map as c, src.image_acc_from_uni as i",
			    "WHERE c.gene_id=i.accession;", sep=" ", collapse="")
		dbGetQuery(db, sql)
	} else if (baseMapType=='ug') {
                message(cat("baseMapType is ug"))
		dbGetQuery(db, "INSERT INTO probe2acc SELECT probe_id, NULL FROM probes_ori;")
		sql <- paste("INSERT INTO probe2gene",
			    "SELECT c.probe_id, u.gene_id",
			    "FROM curr_map as c, src.unigene as u",
			    "WHERE c.gene_id=u.unigene_id;", sep=" ", collapse="")
		dbGetQuery(db, sql)
	} else if (baseMapType=='refseq') {
                message(cat("baseMapType is refseq"))
		dbGetQuery(db, "CREATE INDEX cm1 ON curr_map(probe_id);")
                dbGetQuery(db, "CREATE INDEX cm2 ON curr_map(gene_id);")
                ##just keep ALL of the accessions that went in from the base map (whether or NOT they successfully map over)
		dbGetQuery(db, "INSERT INTO probe2acc SELECT DISTINCT probe_id, gene_id FROM curr_map;")
		sql <- paste("INSERT INTO probe2gene", 
			    "SELECT c.probe_id, a.gene_id",
			    "FROM curr_map as c, src.refseq as a",
			    "WHERE c.gene_id=a.accession;", sep=" ", collapse="") 
		dbGetQuery(db, sql)
	} else 	{ ### type='gb' or 'gbNRef'
                message(cat("baseMapType is gb or gbNRef"))
		dbGetQuery(db, "CREATE INDEX cm1 ON curr_map(probe_id);")
		dbGetQuery(db, "CREATE INDEX cm2 ON curr_map(gene_id);")
                ##just keep ALL of the accessions that went in from the base map (whether or NOT they successfully map over)
		dbGetQuery(db, "INSERT INTO probe2acc SELECT DISTINCT probe_id, gene_id FROM curr_map;")
		sql <- paste("INSERT INTO probe2gene", 
			    "SELECT c.probe_id, a.gene_id",
			    "FROM curr_map as c, src.accession as a",
			    "WHERE c.gene_id=a.accession;", sep=" ", collapse="") 
		dbGetQuery(db, sql)
		sql <- paste("DELETE FROM curr_map WHERE probe_id IN",
			"(SELECT probe_id FROM probe2gene);", sep=" ", collapse="")
 		dbGetQuery(db, sql)
		sql <- paste("INSERT INTO probe2gene ", 
			    "SELECT c.probe_id, a.gene_id ",
			    "FROM curr_map as c, src.accession_unigene as a",  #The HORRIBLY MISNAMED accession_unigene table contains refseq IDs as well as GenBank IDs but no unigene IDs!
                            "WHERE c.gene_id=a.accession;", sep=" ", collapse="") #This name is not my fault!  NOT MY FAULT! - Marc
		dbGetQuery(db, sql)
	}
        
        ##Code for dealing with other sources of IDs (try to find other EG matches and if so, stick them in!)        
	dbGetQuery(db, "DROP TABLE curr_map;")
	lapply(otherSrc, function(thisOtherName) {
            clnVals <- thisOtherName
            sqlIns <- "INSERT INTO other (probe_id,gene_id) VALUES(?,?);"
            dbBeginTransaction(db)
            rset <- dbSendPreparedQuery(db, sqlIns, clnVals)
            dbClearResult(rset)
            dbCommit(db)
	})
	dbGetQuery(db, "DELETE FROM other WHERE gene_id IN ('NA', '');")

        ##I need to know if there is refseq, unigene or GenBank data for the organism in question.
        tables = dbGetQuery(db, "SELECT name FROM src.sqlite_master;")       
        
        ##The following will insert any unigene IDs into the database AS entrez gene IDs, (by joining with the other table)
        if(length(grep("unigene",tables))>0){
            sql <- "INSERT INTO probe2gene SELECT DISTINCT m.probe_id, u.gene_id FROM other as m INNER JOIN src.unigene as u WHERE m.gene_id=u.unigene_id;"
            dbGetQuery(db, sql)
        }
        ##The following will insert any missing refseq IDs into the database AS entrez gene IDs, (by joining with the other table)
        if(length(grep("refseq",tables))>0){
            sql <- "INSERT INTO probe2gene SELECT DISTINCT m.probe_id, r.gene_id FROM other as m INNER JOIN src.refseq as r WHERE m.gene_id=r.accession;"
            dbGetQuery(db, sql)
        }        
        ##The following will insert any missing GenBank IDs into the database AS entrez gene IDs, (by joining with the other table)
        if(length(grep("accession",tables))>0){
            sql <- "INSERT INTO probe2gene SELECT DISTINCT m.probe_id, gb.gene_id FROM other as m INNER JOIN src.accession as gb WHERE m.gene_id=gb.accession;"
            dbGetQuery(db, sql)
        }
        
        ##The following will insert any missing Entrez Gene IDs into the database AS themselves, 
        if(length(grep("EGList",tables))>0){
            sql <- "INSERT INTO probe2gene SELECT DISTINCT probe_id, gene_id FROM other WHERE gene_id IN (SELECT gene_id FROM src.EGList);"
            dbGetQuery(db, sql)
        }

        dbGetQuery(db, "INSERT INTO probe2gene SELECT probe_id, NULL FROM probes_ori WHERE probe_id NOT IN (SELECT probe_id FROM probe2gene);")
	dbGetQuery(db, "CREATE INDEX p1 ON probe2gene(probe_id);")
	dbGetQuery(db, "INSERT INTO temp_probe_map SELECT DISTINCT p.probe_id, p.gene_id, a.gene_id FROM probe2acc as a LEFT OUTER JOIN probe2gene as p ON a.probe_id=p.probe_id;")
        dbGetQuery(db, "DROP TABLE other;")
        dbGetQuery(db, "DROP TABLE probe2gene;")
        dbGetQuery(db, "DROP TABLE probe2acc;")
        dbGetQuery(db, "DROP TABLE probes_ori;")
        
        ##Now I need to attach the chipsrc and use verify that all the gene_IDs are in fact contained in the genes table of the corresponding chipsrc DB...
        dbGetQuery(db, paste("ATTACH '", chipSrc,"' AS chipSrc;", sep="", collapse=""))
        dbGetQuery(db, paste("UPDATE temp_probe_map SET gene_id = NULL WHERE gene_id NOT IN (SELECT gene_id FROM chipSrc.genes);",sep=""))

        
        ##Now I need to decide which probes have multiple mappings and flag those.
        probeData <- dbGetQuery(db, "SELECT * FROM temp_probe_map;")
        modifiedData <- labelDuplicatedProbes(probeData)
        sqlIns <- "INSERT INTO probe_map (probe_id,gene_id,accession,is_multiple) VALUES (?,?,?,?);"
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sqlIns, modifiedData)
        dbClearResult(rset)
        dbCommit(db)

##         ##Drop that temp table
        dbGetQuery(db, "DROP TABLE temp_probe_map;")

        dbDisconnect(db)
	pkgName
}



getMapForBiocChipPkg <- function(csvFileName, pkgName, chipMapSrc, chipSrc, 
				otherSrc=character(0),
				baseMapType="gbNRef", 
				outputDir=".",
				allowHomologyGene=FALSE) {
    baseMaps = makeBaseMaps(csvFileName=csvFileName, targetDir=outputDir)
    
    #pre-clean the otherSrc's and pass them along as lists of lists:  #TODO: this fails if there are not any otherSrc's
    otherSrcObjs = list()
    if(length(otherSrc)==0){
        otherSrc = otherSrcObjs
    }else{
        for(i in 1:length(otherSrc)){
            otherSrcObjs[[i]] = cleanSrcMap(otherSrc[[i]])
        }
    }
    #The 1st item in baseMaps is always the genbank ID
    if (baseMapType == "eg"){
        baseMap <- as.data.frame(baseMaps[[2]])
    } else {
        baseMap <- as.data.frame(baseMaps[[1]])
        otherSrcObjs[[length(otherSrcObjs) + 1]] <- as.data.frame(baseMaps[[2]])
    }
    
    ##just for debugging (disable the rest of the time)
    if(FALSE){
        printOutBaseMaps(baseMap, pkgName, otherSrcObjs)
    }
    
    probe2gene(baseMap=baseMap, 
               baseMapType=baseMapType, 
               otherSrc=otherSrcObjs,
               chipMapSrc=chipMapSrc,
               chipSrc=chipSrc,
               pkgName=pkgName,
               outputDir=outputDir,
               allowHomologyGene=allowHomologyGene)
}

getMapForOtherChipPkg <- function(filePath,
                                  pkgName,
                                  chipMapSrc,
                                  chipSrc,
                                  otherSrc=character(0),
                                  baseMapType="gbNRef",
                                  outputDir=".",
                                  allowHomologyGene=FALSE) {
    #pre-clean the otherSrc's and pass them along as lists of lists:  TODO: fix no otherSrc bug here too
    otherSrcObjs = list()
    if(length(otherSrc)==0){
        otherSrc = otherSrcObjs
    }else{
        for(i in 1:length(otherSrc)){
            otherSrcObjs[[i]] = cleanSrcMap(otherSrc[[i]])
        }
    }
    baseMap = cleanSrcMap(filePath)

    ##just for debugging (disable the rest of the time)
    if(FALSE){
        printOutBaseMaps(baseMap, pkgName, otherSrcObjs)
    }

    probe2gene(baseMap=baseMap,
               baseMapType=baseMapType,
               otherSrc=otherSrcObjs,
               chipMapSrc=chipMapSrc,
               chipSrc=chipSrc,
               pkgName=pkgName,
               outputDir=outputDir,
               allowHomologyGene=allowHomologyGene)
}


getMapForYeastChipPkg <- function(affy, fileName, pkgName, outputDir=".") {

    baseMaps = data.frame()
    if(affy==TRUE){
        baseMaps = makeBaseMaps (csvFileName=fileName,
                        targetDir=outputDir) ##,
                        ##outputPrefix=pkgName)
    }
        require("RSQLite")
	drv <- dbDriver("SQLite")
	outputFile <- file.path(outputDir, paste(pkgName, "sqlite", sep="."))	
	db <- dbConnect(drv, outputFile)
	dbGetQuery(db,
                "CREATE TABLE metadata (name VARCHAR(80) PRIMARY KEY, value VARCHAR(255) );")
        metadata_sql <- paste("INSERT INTO metadata VALUES ('PKGNAME', '",
                                pkgName, "');", sep="", collapse="")
	dbGetQuery(db, metadata_sql)
	dbGetQuery(db, 
		"CREATE TABLE probe_map (probe_id TEXT NOT NULL, systematic_name TEXT, is_multiple SMALLINT NOT NULL);");
    if(affy==TRUE){
        probeData <- as.data.frame(baseMaps[[1]])
        clnVals <- labelDuplicatedProbes(probeData)
        sqlIns <- "INSERT INTO probe_map VALUES(?,?,?);"
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sqlIns, clnVals)
        dbClearResult(rset)
        dbCommit(db)        
    }else
    {
        probeData <- cleanSrcMap(fileName)
        clnVals <- labelDuplicatedProbes(probeData)
        sqlIns <- "INSERT INTO probe_map VALUES(?,?,?);"
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sqlIns, clnVals)
        dbClearResult(rset)
        dbCommit(db)
    }
	dbGetQuery(db,
		"UPDATE probe_map SET systematic_name=NULL WHERE systematic_name='NA';")
	dbDisconnect(db)
    
	pkgName
}


#Need to add fileName and affy param here
getMapForArabidopsisChipPkg <- function(affy, fileName, pkgName, chipMapSrc, outputDir=".") {

    #Note: this function and the associated chipmapsrc database for Arabidopsis both assume that affy will only
    #ever make two arrays for arabidposis.
    #if this changes, then the database and this function will need to be updated (this is unlikely however)
    #for non-affy functions, the map will be read in "as is" and used to make a package.
        require("RSQLite")
	drv <- dbDriver("SQLite")
	outputFile <- file.path(outputDir, paste(pkgName, "sqlite", sep="."))	
	db <- dbConnect(drv, outputFile)
	dbGetQuery(db,
                "CREATE TABLE metadata (name VARCHAR(80) PRIMARY KEY, value VARCHAR(255) );")
        metadata_sql <- paste("INSERT INTO metadata VALUES ('PKGNAME', '",
                                pkgName, "');", sep="", collapse="")
	dbGetQuery(db, metadata_sql)
	dbGetQuery(db, 
		"CREATE TABLE probe_map (probe_id TEXT NOT NULL, gene_id TEXT);");
        dbGetQuery(db, paste("ATTACH DATABASE '",chipMapSrc,"' AS src;",sep=""))

    if(affy==TRUE){#Warning: This will only work if the arabidopsis package is one of "the" TWO affy packages.

	if (pkgName == "ag"){url_name="TAIRAGURL"}
	else if (pkgName == "ath1121501"){url_name="TAIRATHURL"}

        insert_sql <- paste("INSERT INTO metadata SELECT 'TAIRCHIPMAPURL', value FROM src.metadata WHERE name='", url_name, "';", sep="")
	dbGetQuery(db, insert_sql);
                
      	insert_sql <- paste("INSERT INTO probe_map SELECT * FROM src.", pkgName, ";", sep="")
       	dbGetQuery(db, insert_sql);
    }else
    {        
        clnVals <- cleanSrcMap(fileName)
        sqlIns <- "INSERT INTO probe_map VALUES(?,?);"
        dbBeginTransaction(db)
        rset <- dbSendPreparedQuery(db, sqlIns, clnVals)
        dbClearResult(rset)
        dbCommit(db)
        
    }
	dbGetQuery(db, "DETACH src;");			
	dbDisconnect(db)
	pkgName
}


makeUniversalMapping <- function(pkgName,
                                chipSrc,
                                baseMapType="eg",
                                outputDir=".") {

	## The rest of this will just make a map of ALL the EGs
        require("RSQLite")
        drv <- dbDriver("SQLite")
        outputFile <- file.path(outputDir, paste(pkgName, "sqlite", sep="."))
        db <- dbConnect(drv, outputFile)
        dbGetQuery(db,
                "CREATE TABLE metadata (name VARCHAR(80) PRIMARY KEY, value VARCHAR(255) );")
        metadata_sql <- paste("INSERT INTO metadata VALUES ('PKGNAME', '",
                                pkgName, "');", sep="", collapse="")
        dbGetQuery(db, metadata_sql)

	#no source file so we have to grab this from the chipmapsrc
        create_sql <- "CREATE TABLE probe_map (gene_id TEXT NOT NULL);"
        dbGetQuery(db, create_sql)

        attach_sql <- paste("ATTACH DATABASE '",chipSrc,"' AS src;", sep="", collapse="")
        dbGetQuery(db, attach_sql)

        #just grab ALL Entrez genes (without any constraint other than that they be UNIQUE) from chipsrc_XXX
        insert_sql <- "INSERT INTO probe_map SELECT DISTINCT gene_id from src.genes ;"
        dbGetQuery(db, insert_sql)
        
        dbDisconnect(db)
        
        pkgName

}
