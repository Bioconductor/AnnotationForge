\name{@ANNOBJPREFIX@GO}
\alias{@ANNOBJPREFIX@GO}
\alias{@ANNOBJPREFIX@GO2EG}
\alias{@ANNOBJPREFIX@GO2ALLEGS}
\title{Maps between Entrez Gene IDs and Gene Ontology (GO) IDs}
\description{
  @ANNOBJPREFIX@GO is an R object that provides mappings between
  entrez gene identifiers and the GO identifiers that they are 
  directly associated with. This mapping and its reverse mapping do NOT
  associate the child terms from the GO ontology with the gene.  Only
  the directly evidenced terms are represented here.

  @ANNOBJPREFIX@GO2ALLEGS is an R object that provides mappings between
  a given GO identifier and all of the Entrez Gene identifiers annotated
  at that GO term OR TO ONE OF IT'S CHILD NODES in the GO ontology.
  Thus, this mapping is much larger and more inclusive than
  @ANNOBJPREFIX@GO2EG.
}
\details{
  If @ANNOBJPREFIX@GO is cast as a list, each Entrez Gene identifier is
  mapped to a list of lists.  The names on the outer list are GO
  identifiers.  Each inner list consists of three named elements: GOID,
  Ontology, and Evidence.

  The GOID element matches the GO identifier named in the outer list and
  is included for convenience when processing the data using 'lapply'.

  The Ontology element indicates which of the three Gene Ontology
  categories this identifier belongs to.  The categories are biological
  process (BP), cellular component (CC), and molecular function (MF).
     
  The Evidence element contains a code indicating what kind of evidence
  supports the association of the GO identifier to the Entrez Gene
  id. Some of the evidence codes in use include:

  IMP: inferred from mutant phenotype  

  IGI: inferred from genetic interaction

  IPI: inferred from physical interaction  

  ISS: inferred from sequence similarity 

  IDA: inferred from direct assay  

  IEP: inferred from expression pattern  

  IEA: inferred from electronic annotation  

  TAS: traceable author statement  

  NAS: non-traceable author statement  

  ND: no biological data available  

  IC: inferred by curator
  
  
  A more complete listing of evidence codes can be found at:

  \url{http://www.geneontology.org/GO.evidence.shtml}
  
    
  If @ANNOBJPREFIX@GO2ALLEGS or @ANNOBJPREFIX@GO2EG is cast as a list,
  each GO term maps to a named vector of entrez gene identifiers and
  evidence codes. A GO identifier may be mapped to the same entrez gene
  identifier more than once but the evidence code can be
  different. Mappings between Gene Ontology identifiers and Gene
  Ontology terms and other information are available in a separate data
  package named GO.
  
  Whenever any of these mappings are cast as a data.frame, all the
  results will be output in an appropriate tabular form.
  
  Mappings between entrez gene identifiers and GO information were
  obtained through their mappings to Entrez Gene identifiers. NAs are
  assigned to entrez gene identifiers that can not be mapped to any Gene
  Ontology information. Mappings between Gene Ontology identifiers an
  Gene Ontology terms and other information are available in a separate
  data package named GO.
  
  All mappings were based on data provided by: @GOSOURCE@
  
  For GO2ALL style mappings, the intention is to return all the genes
  that are the same kind of term as the parent term (based on the idea
  that they are more specific descriptions of the general term).
  However because of this intent, not all relationship types will be
  counted as offspring for this mapping.  Only "is a" and "has a" style
  relationships indicate that the genes from the child terms would be
  the same kind of thing. 
}
\seealso{
  \code{\link{@ANNOBJPREFIX@GO2ALLEGS}}.
}
\references{
\url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/}
}
\examples{
    x <- @ANNOBJPREFIX@GO
    # Get the entrez gene identifiers that are mapped to a GO ID
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    if(length(xx) > 0) {
        # Try the first one
        got <- xx[[1]]           
        got[[1]][["GOID"]]
        got[[1]][["Ontology"]]
        got[[1]][["Evidence"]]
    }
    # For the reverse map:
    # Convert to a list
    xx <- as.list(@ANNOBJPREFIX@GO2EG)
    if(length(xx) > 0){
        # Gets the entrez gene ids for the top 2nd and 3nd GO identifiers
        goids <- xx[2:3]
        # Gets the entrez gene ids for the first element of goids
        goids[[1]]
        # Evidence code for the mappings
        names(goids[[1]])
    }
    # For @ANNOBJPREFIX@GO2ALLEGS
    xx <- as.list(@ANNOBJPREFIX@GO2ALLEGS)
    if(length(xx) > 0){
        # Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
        goids <- xx[2:3]
        # Gets all the Entrez Gene identifiers for the first element of goids
        goids[[1]]
        # Evidence code for the mappings
        names(goids[[1]])
    }
}
\keyword{datasets}
