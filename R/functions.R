.onLoad <- function(libname, pkgname) {  

  #Load all background information (independent of analysis specific parameters, such as abstraction layer, set of interest, etc.)
  data(mesh_background_frequencies_MEDLINE, mesh_preferred_terms, mesh_preferred_terms_dict, mesh_hierarchy_TC, cofhash110,  background_total, package=pkgname, envir=parent.env(environment()))
  }


#' Abstraction Layer on annotation terms (MeSH)
#' 
#' @description This function creates an Abstraction Layer object that can be used to represent input data (e.g., adverse event terms) and a background set (all MEDLINE) at a given abstraction level.
#'
#' @param ATC_annotations The ATC code for the drug of interest.
#' @param MeSH_annotations The MeSH term for the adverse event of interest.
#' @param range The abstraction layer range based on Information Content of terms to be included.
#' @param fixed The abstraction level is defined by a fixed level in the MeSH tree hierarchcy (e.g. level 2).
#' 
#' @details This function uses hierarchical information from MeSH and Information Content (IC) that has been precalculated for every MeSH term based on aggregated counts on the entire MEDLINE background set to aggregate individual specific terms into more general terms from the specified abstraction level. Alternatively, MeSH terms can be mapped to parent terms at a specific fixed level in the MeSH tree hierarchy (e.g., 2nd level).
#'
#' @return An object containing the abstraction level information for the specified level (list containing the range, the layer, the drugs used as the input parameter, the terms in the annotation file that are covered at the specified abstraction layer, and a set of higher-level terms that represent these terms).
#'
#' @examples {
#'
#'   library("GEA")
#'   
#'   drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt")
#'   events <- AdverseEvents("manifestations_MEDLINE_360k.txt")
#'   abs710    <- getAbstractionLayer(drugs,events,c(7,10))
#' }
#'
#' @export
getAbstractionLayer <- function(ATC_annotations,MeSH_annotations,range=c(),fix=c()) {

  getTargetTermSet <- function(abstraction_range=c(), fix=c(), trees=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")) {
    
    if (length(abstraction_range)>0) {
      
      return(row.names(mesh_background_frequencies_MEDLINE[mesh_background_frequencies_MEDLINE$IC>=abstraction_range[1] & mesh_background_frequencies_MEDLINE$IC<abstraction_range[2]  ,] ))
    }
    
    if (length(fix)>0) {

      #return(as.vector(unique(session$mesh_levels[level==2,"ID" , with=FALSE]  ))$ID)
      return(as.vector(unique(mesh_levels[level==2,"ID" , with=FALSE]  ))$ID)
    }
  } 

  l <- as.data.table(mesh_hierarchy_TC[mesh_hierarchy_TC$parent %in% getTargetTermSet( range, fix),])

  setnames(l, "child", "descriptor")
  setkey(l)

  annotation <- as.data.table(merge(ATC_annotations,MeSH_annotations,by="PMID", allow.cartesian=TRUE))

  annotation[,code:=NULL]
  annotation[,event:=NULL]
  setkey(annotation)

  annotation_anc <-  merge(annotation, l, by  = "descriptor", allow.cartesian=TRUE )  
  setkey(annotation_anc)
  annotation_anc <- unique(annotation_anc)

  alayer <- list("range"=range, "layer"= l, "drugs"=ATC_annotations, "annotation"=unique(annotation), "annotation_anc"=annotation_anc)
  
  return(alayer)
}


#' This function calculates the p-value for GEA.
#'
#' @description This function will detect signals for a given adverse drug event pair based on the biomedical literature (MEDLINE) using GEA.
#'
#' @param atccode The ATC code for the drug of interest.
#' @param event The MeSH term for the adverse event of interest.
#' @param abstraction The abstraction layer object for specified abstraction level.
#' 
#' @details Calculates significane of co-occurrence of drug and adverse event in comparison to large reference set (all MEDLINE). Returns p-value for given drug - adverse event pair.
#'
#' @return P-value for co-occurrence of a given drug - adverse event pair.
#'
#' @examples {
#'
#'   library("GEA")
#'   
#'   drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt")
#'   events <- AdverseEvents("manifestations_MEDLINE_360k.txt")
#'   abs710    <- getAbstractionLayer(drugs,events,c(7,10))
#'   
#'   getPValue("R03DC01","D058186",abs710)
#' }
#'
#' @export
getPValue <- function(drug,event,abstraction) { 

  drugs_i <- abstraction$drugs[code==drug,]
  m <- getEnrichment(drugs_i, abstraction, filter=as.character(event))$annotation_profile
  
  x <- m[m$parent==as.character(event),]$pvalue
  if(length(x)>0) {
    return(x)
  } else {
    return(1)
  }
}  


#' This function calculates the Proportional Reporting Ratio (PRR) for a drug and an adverse event.
#'
#' @description This function will detect signals for given adverse drug event pair based on the biomedical literature (MEDLINE).
#'
#' @param atccode The ATC code for the drug of interest.
#' @param event The MeSH term for the adverse event of interest.
#' @param abstraction The abstraction layer object for specified abstraction level.
#' 
#' @details Constructs a contingency table for given drug and adverse event based on theit occurrrences and co-occurrences in the literature set (provided by the abstraction layer object), aggregates terms to specified abstraction level. Performs zero-cell correction. Returns PRR for given drug - adverse event pair.
#'
#' @return PRR for given drug - adverse event pair
#'
#' @examples {
#'
#'   library("GEA")
#'   
#'   drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt")
#'   events <- AdverseEvents("manifestations_MEDLINE_360k.txt")
#'   abs710    <- getAbstractionLayer(drugs,events,c(7,10))
#'   
#'   getPRR("R03DC01","D058186",abs710)
#' }
#'
#' @export
getPRR <- function(atccode, event, abstraction=getAbstractionLayer(c(4,7))) {
  
  ingred <- abstraction$drugs[ code== atccode ,] 
  ingred <-  unique(ingred$name)[1] 
  
  ab <- nrow(unique(abstraction$annotation_anc[name==ingred,"PMID", with = FALSE]))
  ac <- nrow(unique(abstraction$annotation_anc[parent==event,"PMID", with = FALSE]))
  
  a <- nrow(unique(abstraction$annotation_anc[name==ingred &  parent==event,"PMID",with = FALSE]))
  b<-(ab-a)
  c<-(ac-a)
  d <- nrow(unique(abstraction$annotation_anc[,"PMID", with = FALSE]))-a-b-c
  
  if (a==0 |b==0 | c==0) {
    a=a+.5
    b=b+.5
    c=c+.5
    d=d+.5   
  }
  
  ((a/(a+b))/(c/(c+d)))
}



#' This function executes the enrichment analysis
#'
#' @description This function will analyze which terms are enriched for a given set of interest at the given abstraction level. A list of all terms in the set will be returned together with counts and p-values for enrichment.
#'
#' @param set_of_interest The set of articles for a given set of interest (e.g., a drug, several drugs, a drug class).
#' @param abstraction The abstraction layer object for specified abstraction level.
#' @param threshold Only terms with p-values below threshold are returned.
#' @param filter  Instead for all terms in the set of interest only for the specified terms enrichment will be calculated.
#' @param adjust  If set to TRUE, conditional pvalues will be calculated that adjust for known dependencies among terms.
#' 
#' @details Intersects set of interest with MeSH annotation terms, aggregates terms to specified abstraction level, calculates enrichment for these terms against the general background set (all MEDLINE). All terms that are enrichend (below a given threshold, default = 0.1) will be returned. If filter is set, from all enriched terms only those in the filter set are returned. If the adjust flag is set, enriched terms will be checked for known dependencies based on pairwise conditional hypergeometric test among all enriched terms. This might take considerably longer than the unadjusted EA. 
#'
#' @return Table of enriched terms for given set of interest ordered by p-value.
#'
#' @examples {
#'
#'   library("GEA")
#'   
#'   drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt")
#'   events <- AdverseEvents("manifestations_MEDLINE_360k.txt")
#'   abs710    <- getAbstractionLayer(drugs,events,c(7,10))
#'   
#'   getEnrichment(drugs[ code== "C01EB03" ,] , abs710)$annotation_profile
#' }
#'
#' @export
getEnrichment<- function(set_of_interest, abstraction=getAbstractionLayer(c(4,7)), threshold=0.1, filter=c(),adjust=FALSE) {

  total <- background_total

  abstraction_range <- abstraction$range

  getTerms <- function(ID) {

    cf <- mesh_preferred_terms.dict[[ID]]
    
    if (length(cf)==0) {
      cf = "NA"
    }
    as.vector(cf)
  }

  getFrequency <- function(term) {
    
    mesh_background_frequencies_MEDLINE[term,"count"]
  }

  getPValue <- function(descr,samplesize, observed) {
    
    a <- getFrequency(descr) 
    phyper(observed-1, a , total-a , samplesize, lower.tail=FALSE) 
    
  }

  getAdjustedPValue <- function(l,o,m,p,q,r,z,a1) {
    
    phyper((m-1), l , (a1-l), o, lower.tail=FALSE) * phyper( (r-1), p, (z-p) , q, lower.tail=FALSE)
  }

  
  # create annotation subset: from all papers select those with the entities (e.g., drugs of interest)
  annotation <- abstraction$annotation[PMID %in% set_of_interest$PMID,  ]

  # map to mesh desriptors in the set of interest to parent terms of given abstractionlayer
  annotation_anc1 <- merge(annotation,abstraction$layer,by ="descriptor" ,allow.cartesian=TRUE) 

  # Convert into data.table that contains unique triples of pmid, drugs, and aggregated mesh terms
  annotation_anc  = data.table::copy(annotation_anc1)
  annotation_anc[,descriptor:=NULL]
  annotation_anc[,name:=NULL]
  annotation_anc<-unique(annotation_anc)
  setkey(annotation_anc)
  
  # collect aggregated mesh terms
  annotation_terms <- annotation_anc[, "parent", with = FALSE]
  setnames(annotation_terms, "parent", "parent2")

  # if filter was given as parameter, inersect annotation temrs so far
  # with the ones from the filter
  if(length(filter)>0) {
    
    annotation_terms <- as.data.table(filter)
    setnames(annotation_terms, "filter", "parent2")
    annotation_terms<- annotation_terms[parent2 %in%  annotation_anc$parent ,]
  }
  
  # For each term in set, prepare data needed for enricjment analysis
  # - aggregated term name
  # - total count of papers in this set
  # - count of papers in the set annotated for the aggregated term
  
  parent_term_counts <- unique(as.data.table(annotation_terms))
  parent_term_counts[, total:= nrow(unique(annotation_anc[, "PMID", with = FALSE]))]
  annotation_anc <- annotation_anc[,.(Count=.N), by=parent]
  parent_term_counts[, count:= annotation_anc[parent==parent2, "Count" , with = FALSE]]

  
  # if there is at least one enriched term do enrichment analysis
  # otherwise return empty result set
  if (nrow(parent_term_counts)>0) {

    ##################################################################################
    ### Here calculate un-adjusted p-values for each aggregated term
    ##################################################################################
    parent_term_counts[, pvalue:= getPValue(parent2,total,count)]
  
    # add names for aggregated terms
    parent_term_counts[, term:= lapply( as.vector(parent2) , getTerms)]
    
    # sort by pvalue
    setorder(parent_term_counts,pvalue,-count)

    # Only store terms with enrichmnet above given threshold
    result_matrix <- parent_term_counts[ pvalue<threshold,]
    
    ##################################################
    # Here check if conditional
    ##################################################
    
    #if (nrow(result_matrix)>0) {
    #  result_matrix <- parent_term_counts[pvalue<threshold,,drop=FALSE]
    #}

    
    abstraction <- c( abstraction, list("result_matrix"=result_matrix,
                                        "annotation_profile" = parent_term_counts[order(pvalue,decreasing=FALSE),],
                                        "annotations"= annotation_anc,
                                        "maxp"=matrix(, nrow = 0, ncol = 0),
                                        "dependent"=data.frame(matrix(, nrow = 0, ncol = 0)),
                                        "dependencies"=matrix(, nrow = 0, ncol = 0)))

    return(abstraction)

  } else {
    
        return(list("result_matrix"=data.frame(matrix(, nrow = 0, ncol = 0)),
                   "annotation_profile" = data.frame(matrix(, nrow = 0, ncol = 0)),
                   "annotations"= annotation_anc,
                   "maxp"=matrix(, nrow = 0, ncol = 0),
                   "dependent"=data.frame(matrix(, nrow = 0, ncol = 0)),
                   "dependencies"=matrix(, nrow = 0, ncol = 0)))
  }
}
