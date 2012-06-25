## Function to get all these populate functions out of the namespace.
populateDB = function(schema, ...){
    switch(schema,           
           "HUMANCHIP_DB" = return(AnnotationForge:::popHUMANCHIPDB(...)),
           "MOUSECHIP_DB"  = return(AnnotationForge:::popMOUSECHIPDB(...)),
           "RATCHIP_DB"  = return(AnnotationForge:::popRATCHIPDB(...)),
           "FLYCHIP_DB"  = return(AnnotationForge:::popFLYCHIPDB(...)),
           "YEASTCHIP_DB"  = return(AnnotationForge:::popYEASTCHIPDB(...)),
           "ZEBRAFISHCHIP_DB"  = return(AnnotationForge:::popZEBRAFISHCHIPDB(...)),
           "ECOLICHIP_DB"  = return(AnnotationForge:::popECOLICHIPDB(...)),
           "CANINECHIP_DB"  = return(AnnotationForge:::popCANINECHIPDB(...)),
           "BOVINECHIP_DB"  = return(AnnotationForge:::popBOVINECHIPDB(...)),
           "WORMCHIP_DB"  = return(AnnotationForge:::popWORMCHIPDB(...)),
           "PIGCHIP_DB"  = return(AnnotationForge:::popPIGCHIPDB(...)),
           "CHICKENCHIP_DB"  = return(AnnotationForge:::popCHICKENCHIPDB(...)),
           "XENOPUSCHIP_DB"  = return(AnnotationForge:::popXENOPUSCHIPDB(...)),
           "ARABIDOPSISCHIP_DB"  = return(AnnotationForge:::popARABIDOPSISCHIPDB(...)),
           
           "HUMAN_DB"  = return(AnnotationForge:::popHUMANDB(...)),
           "MALARIA_DB"  = return(AnnotationForge:::popMALARIADB(...)),
           "MOUSE_DB"  = return(AnnotationForge:::popMOUSEDB(...)),
           "RAT_DB"  = return(AnnotationForge:::popRATDB(...)),
           "FLY_DB"  = return(AnnotationForge:::popFLYDB(...)),
           "YEAST_DB"  = return(AnnotationForge:::popYEASTDB(...)),
           "ZEBRAFISH_DB"  = return(AnnotationForge:::popZEBRAFISHDB(...)),
           "CANINE_DB"  = return(AnnotationForge:::popCANINEDB(...)),
           "BOVINE_DB"  = return(AnnotationForge:::popBOVINEDB(...)),
           "WORM_DB"  = return(AnnotationForge:::popWORMDB(...)),
           "PIG_DB"  = return(AnnotationForge:::popPIGDB(...)),
           "CHICKEN_DB"  = return(AnnotationForge:::popCHICKENDB(...)),
           "ARABIDOPSIS_DB"  = return(AnnotationForge:::popARABIDOPSISDB(...)),
           "ECOLI_DB"  = return(AnnotationForge:::popECOLIDB(...)),           
           "CHIMP_DB"  = return(AnnotationForge:::popCHIMPDB(...)),
           "RHESUS_DB"  = return(AnnotationForge:::popRHESUSDB(...)),
           "RHESUSCHIP_DB" = return(AnnotationForge:::popRHESUSCHIPDB(...)),
           "ANOPHELES_DB"  = return(AnnotationForge:::popANOPHELESDB(...)),
           "XENOPUS_DB"  = return(AnnotationForge:::popXENOPUSDB(...))
          )
}


## Function to get all the make***CHIP_DB functions out of the namespace.
makeDBPackage = function(schema, ...){
    switch(schema,           
           "HUMANCHIP_DB" = return(AnnotationForge:::.makeHUMANCHIP_DB(...)),
           "MOUSECHIP_DB"  = return(AnnotationForge:::.makeMOUSECHIP_DB(...)),
           "RATCHIP_DB"  = return(AnnotationForge:::.makeRATCHIP_DB(...)),
           "FLYCHIP_DB"  = return(AnnotationForge:::.makeFLYCHIP_DB(...)),
           "YEASTCHIP_DB"  = return(AnnotationForge:::.makeYEASTCHIP_DB(...)),
           "ZEBRAFISHCHIP_DB"  = return(AnnotationForge:::.makeZEBRAFISHCHIP_DB(...)),
           "ECOLICHIP_DB"  = return(AnnotationForge:::.makeECOLICHIP_DB(...)),
           "CANINECHIP_DB"  = return(AnnotationForge:::.makeCANINECHIP_DB(...)),
           "BOVINECHIP_DB"  = return(AnnotationForge:::.makeBOVINECHIP_DB(...)),
           "WORMCHIP_DB"  = return(AnnotationForge:::.makeWORMCHIP_DB(...)),
           "PIGCHIP_DB"  = return(AnnotationForge:::.makePIGCHIP_DB(...)),
           "CHICKENCHIP_DB"  = return(AnnotationForge:::.makeCHICKENCHIP_DB(...)),
           "XENOPUSCHIP_DB"  = return(AnnotationForge:::.makeXENOPUSCHIP_DB(...)),
           "RHESUSCHIP_DB" = return(AnnotationForge:::.makeRHESUSCHIP_DB(...)),
           "ARABIDOPSISCHIP_DB"  = return(AnnotationForge:::.makeARABIDOPSISCHIP_DB(...))           
          )
}


##Makes a simpleBimap for tables that are added outside of standard AnnotationForge Schemas.
##This function requires that the bimap map from a single table in the DB.
createSimpleBimap <- function(tablename, Lcolname, Rcolname,
                              datacache,
                              objName=as.character(NA),
                              objTarget=as.character(NA))
{
    seed <- list(
                 objName=objName,
                 objTarget=objTarget,
                 Class="AnnDbBimap",
                 L2Rchain=list(
                   list(
                        tablename=tablename,
                        Lcolname=Lcolname,
                        Rcolname=Rcolname
                        )
                   ),
                 datacache=datacache
                 )
    AnnotationDbi:::createAnnDbBimap(seed, list())
}

