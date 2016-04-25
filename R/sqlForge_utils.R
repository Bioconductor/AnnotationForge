## Function to get all these populate functions out of the namespace.
populateDB = function(schema, ...){
    switch(schema, 
           ## ChipDb packages 
           "HUMANCHIP_DB" = return(popHUMANCHIPDB(...)),
           "MOUSECHIP_DB"  = return(popMOUSECHIPDB(...)),
           "RATCHIP_DB"  = return(popRATCHIPDB(...)),
           "FLYCHIP_DB"  = return(popFLYCHIPDB(...)),
           "YEASTCHIP_DB"  = return(popYEASTCHIPDB(...)),
           "ZEBRAFISHCHIP_DB" = 
               return(popZEBRAFISHCHIPDB(...)),
           "ECOLICHIP_DB"  = return(popECOLICHIPDB(...)),
           "CANINECHIP_DB"  = return(popCANINECHIPDB(...)),
           "BOVINECHIP_DB"  = return(popBOVINECHIPDB(...)),
           "WORMCHIP_DB"  = return(popWORMCHIPDB(...)),
           "PIGCHIP_DB"  = return(popPIGCHIPDB(...)),
           "CHICKENCHIP_DB"  = return(popCHICKENCHIPDB(...)),
           "XENOPUSCHIP_DB"  = return(popXENOPUSCHIPDB(...)),
           "ARABIDOPSISCHIP_DB" = return(popARABIDOPSISCHIPDB(...)),

           ## OrgDb packages 
           "HUMAN_DB"  = return(popHUMANDB(...)),
           "MALARIA_DB"  = return(popMALARIADB(...)),
           "MOUSE_DB"  = return(popMOUSEDB(...)),
           "RAT_DB"  = return(popRATDB(...)),
           "FLY_DB"  = return(popFLYDB(...)),
           "YEAST_DB"  = return(popYEASTDB(...)),
           "ZEBRAFISH_DB"  = return(popZEBRAFISHDB(...)),
           "CANINE_DB"  = return(popCANINEDB(...)),
           "BOVINE_DB"  = return(popBOVINEDB(...)),
           "WORM_DB"  = return(popWORMDB(...)),
           "PIG_DB"  = return(popPIGDB(...)),
           "CHICKEN_DB"  = return(popCHICKENDB(...)),
           "ARABIDOPSIS_DB"  = return(popARABIDOPSISDB(...)),
           "ECOLI_DB"  = return(popECOLIDB(...)),
           "CHIMP_DB"  = return(popCHIMPDB(...)),
           "RHESUS_DB"  = return(popRHESUSDB(...)),
           "RHESUSCHIP_DB" = return(popRHESUSCHIPDB(...)),
           "ANOPHELES_DB"  = return(popANOPHELESDB(...)),
           "XENOPUS_DB"  = return(popXENOPUSDB(...))
          )
}

## Function to get all the make***CHIP_DB functions out of the namespace.
makeDBPackage = function(schema, ...){
    switch(schema, 
           "HUMANCHIP_DB" = return(.makeHUMANCHIP_DB(...)),
           "MOUSECHIP_DB"  = return(.makeMOUSECHIP_DB(...)),
           "RATCHIP_DB"  = return(.makeRATCHIP_DB(...)),
           "FLYCHIP_DB"  = return(.makeFLYCHIP_DB(...)),
           "YEASTCHIP_DB"  = return(.makeYEASTCHIP_DB(...)),
           "ZEBRAFISHCHIP_DB"  = return(.makeZEBRAFISHCHIP_DB(...)),
           "ECOLICHIP_DB"  = return(.makeECOLICHIP_DB(...)),
           "CANINECHIP_DB"  = return(.makeCANINECHIP_DB(...)),
           "BOVINECHIP_DB"  = return(.makeBOVINECHIP_DB(...)),
           "WORMCHIP_DB"  = return(.makeWORMCHIP_DB(...)),
           "PIGCHIP_DB"  = return(.makePIGCHIP_DB(...)),
           "CHICKENCHIP_DB"  = return(.makeCHICKENCHIP_DB(...)),
           "XENOPUSCHIP_DB"  = return(.makeXENOPUSCHIP_DB(...)),
           "RHESUSCHIP_DB" = return(.makeRHESUSCHIP_DB(...)),
           "ARABIDOPSISCHIP_DB"  = return(.makeARABIDOPSISCHIP_DB(...))
          )
}
