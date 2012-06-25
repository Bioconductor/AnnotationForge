### R code from vignette source 'SQLForge.Rnw'

###################################################
### code chunk number 1: FileDemo
###################################################
library(AnnotationForge)
read.table(system.file("extdata", "hcg110_ID",
                       package="AnnotationDbi"),
           sep = "\t", header = FALSE, as.is = TRUE)[1:5,]


###################################################
### code chunk number 2: availableDB0s
###################################################
  available.db0pkgs()  


###################################################
### code chunk number 3: GetIntermedDB (eval = FALSE)
###################################################
##   source("http://bioconductor.org/biocLite.R")
##   biocLite("human.db0")


###################################################
### code chunk number 4: PopulateSetup
###################################################
hcg110_IDs = system.file("extdata",
                           "hcg110_ID",
                           package="AnnotationDbi")

myMeta = c("DBSCHEMA"="HUMANCHIP_DB",
           "ORGANISM"="Homo sapiens",
           "SPECIES"="Human",
           "MANUFACTURER"="Affymetrix",
           "CHIPNAME"="Human Cancer G110 Array ",
           "MANUFACTURERURL"="http://www.affymetrix.com")


###################################################
### code chunk number 5: TempDirPopulate
###################################################
tmpout = tempdir()
##To see what chip packages are available:
available.chipdbschemas()
##Then you can make a DB using that schema.
populateDB("HUMANCHIP_DB", affy = FALSE, prefix = "hcg110Test",
    fileName = hcg110_IDs, metaDataSrc = myMeta,
    baseMapType = "gb", outputDir = tmpout)


###################################################
### code chunk number 6: MakeAnnDbPkg
###################################################
seed <- new("AnnDbPkgSeed",             
            Package = "hcg110Test.db",
            Version = "1.0.0",
            PkgTemplate = "HUMANCHIP.DB",
            AnnObjPrefix = "hcg110Test")

makeAnnDbPkg(seed,
             file.path(tmpout, "hcg110Test.sqlite"),
             dest_dir = tmpout)


###################################################
### code chunk number 7: cleanup
###################################################
file.remove(file.path(tmpout, "hcg110Test.sqlite"))
file.rename(file.path(tmpout, "hcg110Test.db"),file.path(tmpout, "fooTest.db"))


###################################################
### code chunk number 8: SQLForge
###################################################
makeDBPackage("HUMANCHIP_DB",
              affy=FALSE,            
              prefix="hcg110",
              fileName=hcg110_IDs,
              baseMapType="gb",
              outputDir = tmpout,
              version="1.0.0",
              manufacturer = "Affymetrix",
              chipName = "Human Cancer G110 Array",
              manufacturerUrl = "http://www.affymetrix.com")


###################################################
### code chunk number 9: cleanup2
###################################################
file.remove(file.path(tmpout, "hcg110.sqlite"))
file.rename(file.path(tmpout, "hcg110.db"),file.path(tmpout, "foo.db"))


###################################################
### code chunk number 10: install (eval = FALSE)
###################################################
## install.packages("packageNameAndPath", repos=NULL, type="source")


###################################################
### code chunk number 11: createSimpleMapping
###################################################
library(hgu95av2.db)
hgu95av2NAMESYMBOL <- createSimpleBimap("gene_info",
                                        "gene_name",
                                        "symbol",
                                        hgu95av2.db:::datacache,
                                        "NAMESYMBOL",
                                        "hgu95av2.db")
##What is the mapping we just made?
hgu95av2NAMESYMBOL
##Display the 1st 4 relationships in this new mapping
as.list(hgu95av2NAMESYMBOL)[1:4]


###################################################
### code chunk number 12: SessionInfo
###################################################
sessionInfo()


