#############################################################################
#############################################################################
###
### AnnDbPkg-maker.R file
###
#############################################################################
#############################################################################


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "AnnDbPkgSeed" class.
###

setClass(
    "AnnDbPkgSeed",
    representation(
        Package="character",            # e.g. "hgu133a2.db"
        Title="character",
        Version="character",            # e.g. "0.0.99"
        License="character", 
        Author="character", 
        Maintainer="character", 
        PkgTemplate="character",        # e.g. "HUMANCHIP.DB"
        DBschema="character",           # e.g. "HUMANCHIP_DB"
        AnnObjPrefix="character",       # e.g. "hgu133a2"
        AnnObjTarget="character",       # e.g. "chip hgu133a2"
        organism="character",           # e.g. "Homo sapiens"
        species="character",            # e.g. "Human"
        manufacturer="character",       # e.g. "Affymetrix"
        chipName="character",           # e.g. "Human Genome U133A 2.0 Array"
        manufacturerUrl="character",    # e.g. "http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133-20"
        biocViews="character"
    ),
    prototype(
        Title=as.character(NA),
        License="Artistic-2.0",
        Author="Marc Carlson",
        Maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
        DBschema=as.character(NA),
        AnnObjPrefix=as.character(NA),
        AnnObjTarget=as.character(NA),
        organism=as.character(NA),
        species=as.character(NA),
        manufacturer=as.character(NA),
        chipName=as.character(NA),
        manufacturerUrl=as.character(NA),
        biocViews=as.character(NA)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some helper functions.
###

initWithDbMetada <- function(x, dbfile)
{
    metadata2slot <- c(
        DBSCHEMA="DBschema",
        ORGANISM="organism",
        SPECIES="species",
        MANUFACTURER="manufacturer",
        CHIPNAME="chipName",
        MANUFACTURERURL="manufacturerUrl"
    )
    dbconn <- dbFileConnect(dbfile)
    on.exit(dbFileDisconnect(dbconn))
    metadata <- dbGetTable(dbconn, "metadata")
    if (!identical(colnames(metadata), c("name", "value")))
        stop("\"metadata\" table has unexpected col names")
    if (any(duplicated(metadata$name))) {
        stop("col \"name\" in \"metadata\" table has duplicated values\n",
             "  (this would never happen if \"name\" was defined as a PRIMARY KEY!)")
    }
    row.names(metadata) <- metadata$name
    for (i in seq_len(length(metadata2slot))) {
        metadata_name <- names(metadata2slot)[i]
        if (!(metadata_name %in% row.names(metadata))) {
            if (metadata_name == "DBSCHEMA")
                stop("'DBSCHEMA' not found in \"metadata\" table")
            next
        }
        slot_name <- metadata2slot[i]
        val <- metadata[metadata_name, "value"]
        if (is.na(slot(x, slot_name))) {
            slot(x, slot_name) <- val
            next
        }
        if (slot(x, slot_name) != val)
            stop(metadata_name, " specified in '", dbfile, "' (\"", val, "\") ",
                 "doesn't match 'x@", slot_name, "' (\"", slot(x, slot_name), "\")")
    }
    if (is.na(x@manufacturerUrl)) {
        x@manufacturerUrl <- ""
        warning("no manufacturerUrl for package ", x@Package)
    }
    x
}

initComputedSlots <- function(x)
{
    if (is.na(x@AnnObjPrefix))
        stop("'AnnObjPrefix' slot must be set for package ", x@Package)
    ## Automatic default for "AnnObjTarget" slot
    if (is.na(x@AnnObjTarget))
        x@AnnObjTarget <- paste("chip", x@AnnObjPrefix)
    ## Automatic default for "Title" slot
    if (is.na(x@Title)) {
        if (is.na(x@manufacturer) || is.na(x@chipName) || is.na(x@AnnObjTarget)) {
            warning("not enough information to set the 'Title' slot for package ", x@Package)
        } else {
            x@Title <- paste(x@manufacturer,
                             " ",
                             x@chipName,
                             " annotation data (",
                             x@AnnObjTarget,
                             ")", sep="")
        }
    } 
    ## Automatic default for "biocViews" slot
    if (is.na(x@biocViews)
     && !is.na(x@organism)
     && !is.na(x@manufacturer)) {
        chip_view <- paste(x@manufacturer, "Chip", sep="")
        org_view <- chartr(" ", "_", x@organism)
        x@biocViews <- paste("AnnotationData", chip_view, org_view,
                             x@AnnObjPrefix, sep=", ")
    }
    x
}

initWithDbDoc <- function(dbfile)
{
    dbconn <- dbFileConnect(dbfile)
    on.exit(dbFileDisconnect(dbconn))
    if(dbExistsTable(dbconn, "map_metadata")){
        map_metadata <- dbGetTable(dbconn, "map_metadata")
        return(map_metadata)
    } else {
        return(NULL)
    }
}

getSymbolValuesForManPages <- function(map_names, dbfile)
{
    map_metadata <- initWithDbDoc(dbfile)
    if(is.null(map_metadata)) return(NULL)
    map_source <- sapply(map_names,
                         function(this_map)
                         {
                             map_index <- which(map_metadata$map_name == this_map)
                             if (length(map_index) > 0) {
                                 this_source <- paste(
                                     map_metadata[map_index, "source_name"],
                                     " \n ",
                                     map_metadata[map_index, "source_url"],
                                     " \n  With a date stamp from the source of:",
                                     map_metadata[map_index, "source_date"],
                                     sep=" ", collapse=" and ")
                             } else {
                                 this_source <- NA
                             }
                             this_source
                         })
    map_source <- gsub("_", "\\_", map_source, fixed=TRUE)
    names(map_source) <- paste(map_names, "SOURCE", sep="")
    map_source
}

removeCommentsFromFile <- function(infile, outfile)
{
    if (!is.character(infile) || length(infile) != 1 || is.na(infile))
        stop("'infile' must be a character string naming a file")
    if (!is.character(outfile) || length(outfile) != 1 || is.na(outfile))
        stop("'outfile' must be a character string naming a file")
    if (file.exists(outfile))
        stop("file '", outfile, "' already exists")
    infile <- file(infile, "r")
    #on.exit(close(infile))
    outfile <- file(outfile, "w")
    #on.exit(close(outfile)) # doesn't seem to work
    while (TRUE) {
        text <- readLines(infile, n=1)
        if (length(text) == 0)
            break
        if (substr(text, 1, 1) != "#")
            writeLines(text, outfile)
    }
    close(outfile)
    close(infile)
}

loadAnnDbPkgIndex <- function(file)
{
    if (missing(file)) {
        file <- system.file("extdata", "GentlemanLab", "ANNDBPKG-INDEX.TXT",
                            package="AnnotationForge")
    } else {
        if (!is.character(file) || length(file) != 1 || is.na(file))
            stop("'file' must be a character string naming a file")
    }
    tmp_file <- file.path(tempdir(), paste(basename(file), "tmp", sep="."))
    removeCommentsFromFile(file, tmp_file)
    index <- read.dcf(tmp_file)
    file.remove(tmp_file)
    index
}





### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers for filtering out innapropriate manual pages from a template.
### 

## This function takes the seed and lists the Mappings
listMappings <- function(x, type){
  ## get seeds
  schema <- x@DBschema ## schema will be like "HUMANCHIP_DB" or "HUMAN_DB"
  if(type=="ChipDb"){
    orgDbName <- getOrgPkgForSchema(schema)
    allSeeds <- NCBICHIP_DB_SeedGenerator(orgDbName)
  }else if(type=="OrgDb"){
    allSeeds <- AnnotationDbi:::NCBIORG_DB_SeedGenerator()
  }
  seeds <- filterSeeds(allSeeds, schema, type)
  ## Then get the names
  unlist(lapply(seeds, function(x){return(x$objName)}))
}


## This function will translate the from the actual mappings to the requisite
## man pages.  The whole point is to get rid of things like flybase from
## humans etc.

## ALSO problematic: CHRLENGTHS.Rd, GO2ALLEGS.Rd UCSCGENES.Rd

filterManPages <- function(doc_template_names, maps, x){
  docs <- sub("\\.Rd$", "", doc_template_names)
  docs <- docs[docs %in% maps]
  ## Add things that will always be needed but are not themselves really bimaps
  docs <- c(docs, "_dbconn" ,"BASE","ORGANISM","MAPCOUNTS")  
  if(!any(c("ECOLI_DB","XENOPUS_DB","ECOLICHIP_DB","XENOPUSCHIP_DB","PIG_DB"
            ,"PIGCHIP_DB") %in% x@DBschema)){
    docs <- c(docs, "CHRLENGTHS")
  }
  paste(docs, ".Rd", sep="")
}


## And I need a wrapper function to help me filter out things that are not in
## the manList when I call createPackage.

.createAnnotPackage <-function(pkgname,destinationDir,originDir,symbolValues,
                              manList, unlink=FALSE, quiet=FALSE){
  
  tdir <- file.path("TEMPANNOTPACKAGEDIRFORFILTERING")
  dir.create(tdir)
#  tdir <- file.path(tempdir()) ## tempdir() causes strange errors...  :(
  file.copy(from = dir(originDir, full.names = TRUE),
            to = tdir,
            recursive = TRUE)
  
  ## Then unlink unwanted man pages from tdir
  manDir <- file.path(tdir, "man")
  manFiles <- dir(manDir)
  rmFiles <- manFiles[!(manFiles %in% manList)]
  rmFiles <- file.path(manDir, rmFiles)
  unlink(rmFiles)
  
  ## Then call createPackage
  createPackage(pkgname=pkgname,
                destinationDir=destinationDir,
                originDir=tdir,
                symbolValues=symbolValues,
                unlink=unlink,
                quiet=quiet)
  ## Then remove our terrible temp dir
  unlink(tdir, recursive=TRUE)
  ## Will need to return to tempdir() if we ever want to be able to do more
  ## than one at a time...  :(
}



## TESTING:
## library(AnnotationForge)
## debug(AnnotationForge:::.createAnnotPackage)
## debug(AnnotationForge:::.makeAnnDbPkg) ## this one is always called.
## debug(AnnotationForge:::.makeAnnDbPkgs) ## This one is called 1st for mine
## debug(AnnotationForge:::.makeAnnDbPkgList) ## called for others
## source("~/proj/Rpacks/AnnotationForge/inst/extdata/GentlemanLab/org-batch-script.R")
## source("~/proj/Rpacks/AnnotationForge/inst/extdata/GentlemanLab/chip-batch-script.R")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "makeAnnDbPkg" new generic.
###

setGeneric("makeAnnDbPkg", signature="x",
    function(x, dbfile, dest_dir=".", no.man=FALSE, ...) 
        standardGeneric("makeAnnDbPkg")
)

## helper to extract metadata
.getOrgDepFromMetadata <- function(dbfile){
    con <- dbConnect(SQLite(), dbfile)
    dbGetQuery(con, "SELECT value FROM metadata WHERE name = 'ORGPKGDEP'")
}

.makeAnnDbPkg <- function(x, dbfile, dest_dir=".", no.man=FALSE, ...){
  x <- initWithDbMetada(x, dbfile)
  x <- initComputedSlots(x)
  dbfile_basename <- basename(dbfile)
  if (dbfile_basename != paste(x@AnnObjPrefix, ".sqlite", sep=""))
    stop("'", dbfile, "': File name doesn't match 'x@AnnObjPrefix' (", x@AnnObjPrefix, ")")
  if (!grepl("^/", x@PkgTemplate)[1]) { ##TODO: this regex seems hacky?
    template_path <- system.file("AnnDbPkg-templates",
                                 x@PkgTemplate,
                                 package="AnnotationForge")
  } else {
    template_path <- x@PkgTemplate
  }
  ann_dbi_version <- installed.packages()['AnnotationDbi','Version']
  ## only define 'org_version' if we are making a chipDb package.
  ## Otherwise it will only cause trouble.
  
  con1 <- dbConnect(dbDriver("SQLite"), dbfile)
  type <- dbGetQuery(con1, 
                     "SELECT value FROM metadata WHERE name='Db type'")
  if(type=="ChipDb"){
    org_version <- installed.packages()['org.Hs.eg.db','Version']
    ## NOCHIPSCHEMA DBs know who they depend on
    if(x@DBschema=="NOCHIPSCHEMA_DB"){
        org_pkg <- as.character(.getOrgDepFromMetadata(dbfile))
    }else{
        org_pkg <- paste0(getOrgPkgForSchema(x@DBschema),".db")
    }
  }else{
    org_version <- "no org version date"
    org_pkg <- "no org pkg required"
  }
  symvals <- list(
                  DBSCHEMA=x@DBschema,
                  PKGTITLE=x@Title,
                  ANNOBJPREFIX=x@AnnObjPrefix,
                  ANNOBJTARGET=x@AnnObjTarget,
                  ORGANISM=x@organism,
                  SPECIES=x@species,
                  MANUF=x@manufacturer,
                  CHIPNAME=x@chipName,
                  MANUFURL=x@manufacturerUrl,
                  AUTHOR=x@Author,
                  MAINTAINER=x@Maintainer,
                  PKGVERSION=x@Version,
                  LIC=x@License,
                  BIOCVIEWS=x@biocViews,
                  DBFILE=dbfile_basename,
                  ANNDBIVERSION=ann_dbi_version,
                  ORGVERSION=org_version,
                  ORGPKGDEP=org_pkg
                  )
  man_dir <- file.path(template_path, "man")
  doc_template_names <- list.files(man_dir, "\\.Rd$")
  if (file.exists(man_dir)) {
    if (!no.man) {
      #is_static <- doc_template_names %in% c("_dbconn.Rd", "_dbfile.Rd")
      #doc_template_names <- doc_template_names[!is_static]


      ## Do this only if your schema is an NCBI* one.
      if(grepl("NCBI",x@PkgTemplate)){
        ## extract the map_names from the bimap definitions
        map_names <- listMappings(x, type)
        ## now use this info to filter to relevant mappings
        doc_template_names <- filterManPages(doc_template_names,
                                             maps=map_names,x)
      }else{## if old school, just use the man pages in template
        map_names <- sub("\\.Rd$", "", doc_template_names)        
      }
      
      if (length(map_names) != 0)
        symvals <- c(symvals, getSymbolValuesForManPages(map_names, dbfile))
    } else {
      doc_template_names <- list()
      unlink(man_dir, recursive=TRUE) # delete template
    }
  }
  if (any(duplicated(names(symvals)))) {
    str(symvals)
    stop("'symvals' contains duplicated symbols (see above)")
  }
  ## Remove NA values
  symvals <- symvals[!sapply(symvals, is.na)]
  .createAnnotPackage(x@Package,
                     destinationDir=dest_dir,
                     originDir=template_path,
                     symbolValues=symvals,
                     manList=doc_template_names)
  ## rename Rd files (prepend the pkg name)
  ## Here is also where we put the man files into the package (after renaming them)
  if (file.exists(man_dir) && !no.man && length(doc_template_names) != 0) {
    doc_path <- file.path(dest_dir, x@Package, "man")
    from_doc_names <- paste(doc_path, doc_template_names, sep=.Platform$file.sep)
    to_doc_names <- paste(x@AnnObjPrefix, doc_template_names, sep="")
    to_doc_names <- paste(doc_path, to_doc_names, sep=.Platform$file.sep)
    mapply(file.rename, from_doc_names, to_doc_names)
  }
  
  dest_db_dir <- file.path(dest_dir, x@Package, "inst", "extdata")
  if (!file.exists(dest_db_dir)
      && !dir.create(dest_db_dir, recursive=TRUE))
    stop("unable to create dest db dir ", dest_db_dir)
  dest_dbfile <- file.path(dest_db_dir, dbfile_basename)
  if (!file.copy(dbfile, dest_dbfile))
    stop("cannot copy file '", dbfile, "' to '", dest_dbfile, "'")
  if(.Platform$OS.type != 'windows'){
    command <- paste("chmod 444", dest_dbfile)
    if (system(command) != 0)
      warning(command, " failed")
  }
  return(invisible(TRUE))
}

setMethod("makeAnnDbPkg", "AnnDbPkgSeed",
          function(x, dbfile, dest_dir=".", no.man=FALSE, ...){
            .makeAnnDbPkg(x, dbfile, dest_dir=dest_dir, no.man=no.man, ...)
          }
          )
        

.makeAnnDbPkgList <- function(x, dbfile, dest_dir=".", no.man=FALSE, ...){
  x$Class <- "AnnDbPkgSeed"
  y <- do.call(new, x)
  makeAnnDbPkg(y, dbfile, dest_dir, no.man)
}

setMethod("makeAnnDbPkg", "list",
    function(x, dbfile, dest_dir=".", no.man=FALSE, ...) {
            .makeAnnDbPkgList(x, dbfile, dest_dir=dest_dir, no.man=no.man, ...)
    }
)

.makeAnnDbPkgs <- function(x, dbfile, dest_dir=".", no.man=FALSE, ...){
  if (missing(dbfile)) {
    dbfile <- system.file("extdata", "GentlemanLab", "ANNDBPKG-INDEX.TXT",
                          package="AnnotationForge")
  }
  index <- loadAnnDbPkgIndex(dbfile)
  if (length(x) != 1) {
    ii <- match(x, index[ , "Package"])
    if (any(is.na(ii)))
      stop("packages ", paste(x[is.na(ii)], collapse=", "),
           " not in ", dbfile)
    index <- index[ii, , drop=FALSE]
  } else if (!is.na(x) && x != "") {
    pkgname <- paste("^", x, "$", sep="")
    ii <- grep(pkgname, index[ , "Package"])
    index <- index[ii, , drop=FALSE]
  }
  filter <- list(...)
  for (j in seq_len(length(filter))) {
    colname <- names(filter)[j]
    if (!(colname %in% colnames(index)))
      stop("unknown field '", colname, "'")
    colvals <- filter[[j]]
    if (!is.character(colvals))
      stop("extra arg values must be of type character")
    index <- index[index[ , colname] %in% colvals, , drop=FALSE]
  }
  pkgnames_in1string <- paste(index[, "Package"], collapse=", ")
  cat(nrow(index), " package(s) to make: ",
      pkgnames_in1string, "\n", sep="")
  for (i in seq_len(nrow(index))) {
    y <- index[i, ]
    y <- as.list(y[!is.na(y)])
    cat("[", i, "/", nrow(index), "] making package ",
        y[["Package"]], ": ", sep="")
    dbfile <- y[["DBfile"]]
    y <- y[names(y) != "DBfile"]
    makeAnnDbPkg(y, dbfile, dest_dir, no.man)
  }
  cat("DONE (", nrow(index), " package(s) made under the ",
      dest_dir, " directory)\n", sep="")
}

setMethod("makeAnnDbPkg", "character",
          function(x, dbfile, dest_dir=".", no.man=FALSE, ...){
            .makeAnnDbPkgs(x, dbfile, dest_dir=dest_dir, no.man=no.man, ...)
          }
          )

