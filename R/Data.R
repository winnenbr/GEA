#' GS
#'
#' @description This function creates an Gold Standard object that can be used to evaluate the performance of drug adverse event signal detection methods.
#'
#' @param gs_name     The name of the gold standard file
#' 
#' @examples {
#' 
#'    gs <- GS("OMOP_GS_drugs.csv")
#' 
#' } 
#' 
#' @return An object containing the Gold Standard data.
#' 
#' @export
GS <- function(gs_name) {
  
  set <- read.csv(system.file("extdata", gs_name, package = "GEA2"))
  
  set$set_bin = NA 
  set <- within(set, set_bin[set=="POS"] <- 1)
  set <- within(set, set_bin[set=="NEG"] <- 0)
  return(as.data.table(set))
  
}


#' Drugs
#'
#' @description This function creates a data table with the drug candidates extracted from MEDLINE
#'
#' @param drug_filename The name of the flat file containing the drug specific information extracted from MEDLINE
#' @param classLevel If TRUE, represent drugs at the drug class level, e.g. "bupivacaine" as "amides" 
#' 
#' @examples {
#' 
#'    drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt")
#'    
#'    #will return the PubMed IDs for the 758 papers in MEDLINE mentioning the drug bupivacaine together with an adverse event
#'    drugs[name=="bupivacaine",]$PMID
#'     
#' } 
#' 
#' @return An object containing for drug candidates: names, ATC codes, and MEDLINE Ids for papers mentioning the drugs potentially in the context of adverse events
#' 
#' @export
Drugs <- function(drug_filename, classLevel=FALSE) {
  
  ATC_annotations <- read.csv(unz(system.file("extdata", gsub(".txt", ".zip", drug_filename), package = "GEA2"), drug_filename ),
    "",
    sep="\t", header=FALSE,check.names=FALSE, fill=TRUE)
  

  ATC_annotations <- ATC_annotations[ATC_annotations$V2=="i",]

  if (classLevel==FALSE) {
    
    ATC_annotations <- ATC_annotations[,c(1,9,10)]
  } else {

    ATC_annotations <- ATC_annotations[,c(1,11,12)]
  }

  colnames(ATC_annotations) <- c("PMID","name","code")
  ATC_annotations <- as.data.table(ATC_annotations[ATC_annotations$code!="",])
  setkey(ATC_annotations,NULL)
  
  return(unique(ATC_annotations))
}


#' AdverseEvents
#'
#' @description This function creates a data table with the event candidates extracted from MEDLINE
#'
#' @param event_filename The name of the flat file containing the event specific information extracted from MEDLINE
#' 
#' @examples {
#' 
#'    events <- AdverseEvents("manifestations_MEDLINE_360k.txt")
#'    
#'    #will return the PubMed IDs for the 933 papers in MEDLINE mentioning Epilepsy as an adverse event
#'    events[event=="Epilepsy",]$PMID
#'     
#' } 
#' 
#' @return An object containing for adverse event candidates: names, MeSH codes, and MEDLINE Ids for papers mentioning potentially adverse events
#' 
#' @export
AdverseEvents <- function(event_filename ) {
  
  MeSH_annotations <- fread(
    system.file("extdata", event_filename, package = "GEA2"),
    sep="\t", header=FALSE, data.table=FALSE
  )
  
  MeSH_annotations <-  MeSH_annotations[,c(1,2,3)]
  colnames(MeSH_annotations) <- c("PMID","event","descriptor" )
  MeSH_annotations<-as.data.table(MeSH_annotations)
  setkey(MeSH_annotations,NULL)
  
  return(unique(MeSH_annotations))
}
