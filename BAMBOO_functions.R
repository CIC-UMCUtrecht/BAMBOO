suppressWarnings(suppressMessages(require(OlinkAnalyze))) #
# suppressWarnings(suppressMessages(require(RUVSeq))) #
# suppressWarnings(suppressMessages(require(ggfortify, verbose = F)))
suppressWarnings(suppressMessages(require(tidyverse))) #
suppressWarnings(suppressMessages(require(pheatmap))) #
suppressWarnings(suppressMessages(require(RColorBrewer))) # 
# suppressWarnings(suppressMessages(require(dbscan))) # 
suppressWarnings(suppressMessages(require(readr))) #
library(preprocessCore) #
suppressWarnings(suppressMessages(require(ggplot2))) #

recode <- dplyr::recode
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
group_by <- dplyr::group_by
case_when <- dplyr::case_when

message("loading the latest function y'all")

BAMBOO_normalization <- function(plateReference, plateSubject, BCs, LODthreshold = 8){
  #plat
  #
  #
  #
  #
  #
  #
  mutate <- dplyr::mutate
  if(is.null(BCs)){
    
    BCs <- intersect(plateReference$SampleID, plateSubject$SampleID)
    
  }
  
  BCs <- removeOutliers(plateReference, plateSubject, BCs, quantileThreshold = 0.99) #based on interquantile outlier detection
  
  flaggedAssays <- flagAssay(plateReference = plateReference, plateSubject = plateSubject, BCs = BCs, BCsAboveLOD = LODthreshold, correlationThreshold = 0)
  
  flaggedCorrelationAssays <- flagAssayCorrelation(plateReference = plateReference, plateSubject = plateSubject, BCs = BCs, BCsAboveLOD = LODthreshold, correlationThreshold = 0)
  
  # flaggedAssays <- NA
  plateReferenceBCs <- plateReference %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays) %>% 
    group_by(Assay) %>% 
    filter(SampleID%in%BCs) %>% # have to decide what to do with the below LOD samples... 
    mutate(pseudoReference = median(NPX, na.rm = T)) %>% #first take the median for the intensity normalized values
    mutate(NPX = case_when(NPX >= LOD ~ NPX)) %>% 
    mutate(BCreference = median(NPX, na.rm = T),
           pseudoReference = case_when(
             AssayFlag ~ pseudoReference,
             T ~ BCreference)
    ) 
  
  plateSubjectBCs <- plateSubject %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays) %>% 
    group_by(Assay) %>% 
    filter(SampleID%in%BCs) %>% # have to decide what to do with the below LOD samples... 
    mutate(pseudoReference = median(NPX, na.rm = T)) %>% #first take the median for the intensity normalized values
    mutate(NPX = case_when(NPX >= LOD ~ NPX)) %>% 
    mutate(BCreference = median(NPX, na.rm = T),
           pseudoReference = case_when(
             AssayFlag ~ pseudoReference,
             T ~ BCreference)
    ) 
  
  mergedPlates <- inner_join(plateReferenceBCs, plateSubjectBCs, by = c("SampleID", "Assay"), suffix = c(".ref", ".sub")) 
  
  if(!all(mergedPlates$SampleID%in%BCs)){
    
    message("sample ", setdiff(mergedPlates$SampleID, BCs), "not present in both plates")
    
  }
  
  modelCoeff <- lmrob(pseudoReference.ref ~ pseudoReference.sub, data = mergedPlates)$coefficients
  
  Adj_factors <- mergedPlates %>% 
    group_by(Assay) %>% 
    mutate(Adj_factors = median(NPX.ref - (NPX.sub*modelCoeff[["pseudoReference.sub"]] + modelCoeff[["(Intercept)"]]), na.rm = T)) %>% 
    mutate(Adj_factors = 
             case_when(
               Assay%in%flaggedAssays ~ pseudoReference.ref - (pseudoReference.sub*modelCoeff[["pseudoReference.sub"]] + modelCoeff[["(Intercept)"]]), 
               T ~ Adj_factors
             )
    ) %>% 
    mutate(Adj_factors = 
             case_when(
               is.na(Adj_factors) & !is.na(pseudoReference.sub) & !is.na(pseudoReference.ref) ~ pseudoReference.ref - (pseudoReference.sub*modelCoeff[["pseudoReference.sub"]] + modelCoeff[["(Intercept)"]]),
               is.na(Adj_factors) ~ LOD.ref - LOD.sub, 
               T ~ Adj_factors)
    ) %>% 
    select(c("Assay", "Adj_factors")) %>% 
    unique()
  
  Adj_factors %>% select(Assay, Adj_factors) %>% distinct()
  
  plateSubject <- plateSubject %>% 
    left_join(Adj_factors, by = "Assay") %>% 
    mutate(NPX = NPX*modelCoeff[["pseudoReference.sub"]] + modelCoeff[["(Intercept)"]] + Adj_factors,
           LOD = LOD*modelCoeff[["pseudoReference.sub"]] + modelCoeff[["(Intercept)"]] + Adj_factors) %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays, AssayFlagCorrelation = Assay%in%flaggedCorrelationAssays)
  
  return(plateSubject)
}