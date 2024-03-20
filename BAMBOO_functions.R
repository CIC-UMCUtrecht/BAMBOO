suppressWarnings(suppressMessages(require(OlinkAnalyze))) #
# suppressWarnings(suppressMessages(require(RUVSeq))) #
# suppressWarnings(suppressMessages(require(ggfortify, verbose = F)))
suppressWarnings(suppressMessages(require(tidyverse))) #
# suppressWarnings(suppressMessages(require(pheatmap))) #
# suppressWarnings(suppressMessages(require(RColorBrewer))) # 
# suppressWarnings(suppressMessages(require(dbscan))) # 
suppressWarnings(suppressMessages(require(readr))) #
# library(preprocessCore) #
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

removeOutliers <- function(plateReference, plateSubject, BCs, quantileThreshold = 1){ #inter Quantile outlier detection -> can be improved
  
  ss <- sumOfSquares(plateReference, plateSubject, samplesOfInterest = BCs, na.rm = T)
  
  outlierSamples <- ss %>% filter(SS > quantile(ss$SS, quantileThreshold)) %>% pull(SampleID)
  
  message(paste0("identiefied the following sample as outlier: ", outlierSamples, "\n"))
  
  return(BCs[!BCs%in%outlierSamples])
  
}

readNPX <- function(file) {
  #function that reads the NPX files using Olink read_NPX function, and gives an error message if it doesn't work
  out <- tryCatch(
    {
      read_NPX(file) 
    },
    error=function(cond) {
      message(paste("File caused a error:", file))
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("File caused a warning:", file))
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}

writeNPX <- function(plateList, plateName, normalizationMethod = "BAMBOO", directory){
  plate <- plateList[[plateName]]
  
  headerAssayInfo <- plate %>% select(c(Panel, Assay, UniProt, OlinkID)) %>% unique() 
  header <- matrix(data = NA, nrow = 2, ncol = nrow(headerAssayInfo) + 5)
  header[1,2] <- "BAMBOO_export_v1"
  header[2,1] <- "NPX data"
  headerAssayInfo1 <- bind_rows(c("Panel" = "Panel", "Assay" = "Assay", "UniProt" = "Uniprot ID", "OlinkID" = "OlinkID"), headerAssayInfo) #%>% mutate(emptyRow = NA) #%>% t()
  headerAssayInfo2 <- data.frame(Panel = unique(plate$Panel),
                                 Assay = c("Plate ID", "QC Warning", "QC Deviation from median", "QC Deviation from median"),
                                 UniProt = c(NA, NA, "Inc Ctrl 2", "Det Ctrl"),
                                 OlinkID = NA)
  
  
  headerAssayInfo <- bind_rows(headerAssayInfo1, headerAssayInfo2) %>% mutate(emptyRow = NA) %>% t() 
  rownames(headerAssayInfo) <- NULL
  header <- rbind(header, headerAssayInfo)
  
  NPXbody <- plate %>% 
    mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% 
    select(c(SampleID, NPX, Assay, PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`)) %>% 
    pivot_wider(names_from = c(Assay), values_from = NPX) %>% 
    relocate(c(PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`), .after = `CSF-1`) %>% as.matrix()
  
  if(!all(colnames(NPXbody)[2:93] == header[4,2:93])){
    message("ABORT ABORT ABORT", plateName)
  }
  
  tail <- plate %>% group_by(Assay) %>% mutate(belowLOD = paste0(mean(NPX < LOD, na.rm = T)*100, "%"), normalization = "BAMBOO") %>% distinct(Assay, LOD, belowLOD, normalization, .keep_all = T) %>% select(Assay, LOD, belowLOD, normalization, AssayFlag) %>% filter(LOD == max(LOD)) %>% mutate(LOD = as.character(LOD)) %>% 
    arrange(match(Assay, header[4,2:93])) %>% t()
  
  if(!all(tail[1,] == header[4,2:93])){
    message("ABORT ABORT ABORT", plateName)
  }
  
  
  
  tail <- cbind(c(NA, "LOD", "Missing Data freq.", "Normalization", "below limit of detection adjustment"), tail)
  
  emptyMat <- matrix(nrow = 5, ncol = 4)
  
  tail <- cbind(tail,emptyMat)
  
  tail[1,] <- tail["AssayFlag",]
  tail["AssayFlag",] <- NA
  
  row.names(tail) <- NULL
  
  
  if(all(is.na(plate$`QC Deviation Det Ctrl`))){
    header <- header[,-c(96,97)]
    NPXbody <- NPXbody[,-c(96,97)]
    tail <- tail[,-c(96,97)]
  }
  
  top <- rbind(header, NPXbody)
  
  exportPlate <- rbind(top, tail)
  colnames(exportPlate) <- NULL
  
  if(!dir.exists(directory)){
    dir.create(directory)
  }
  
  file <- paste0(directory, plateName, "_", normalizationMethod, ".xlsx")
  
  openxlsx::write.xlsx(data.frame(exportPlate), file = file, colNames = F)
  read_NPX(file)
}

writeNPXorganDamage <- function(plateList, plateName, normalizationMethod = "BAMBOO", directory){
  plate <- plateList[[plateName]]
  
  headerAssayInfo <- plate %>% select(c(Panel, Assay, UniProt, OlinkID)) %>% unique() 
  header <- matrix(data = NA, nrow = 2, ncol = nrow(headerAssayInfo) + 5)
  header[1,2] <- "BAMBOO_export_v1"
  header[2,1] <- "NPX data"
  headerAssayInfo1 <- bind_rows(c("Panel" = "Panel", "Assay" = "Assay", "UniProt" = "Uniprot ID", "OlinkID" = "OlinkID"), headerAssayInfo) #%>% mutate(emptyRow = NA) #%>% t()
  headerAssayInfo2 <- data.frame(Panel = unique(plate$Panel),
                                 Assay = c("Plate ID", "QC Warning", "QC Deviation from median", "QC Deviation from median"),
                                 UniProt = c(NA, NA, "Inc Ctrl 2", "Det Ctrl"),
                                 OlinkID = NA)
  
  
  headerAssayInfo <- bind_rows(headerAssayInfo1, headerAssayInfo2) %>% mutate(emptyRow = NA) %>% t() 
  rownames(headerAssayInfo) <- NULL
  header <- rbind(header, headerAssayInfo)
  
  NPXbody <- plate %>% 
    mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% 
    select(c(SampleID, NPX, Assay, PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`)) %>% 
    pivot_wider(names_from = c(Assay), values_from = NPX) %>% 
    relocate(c(PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`), .after = `CALR`) %>% 
    as.matrix()
  
  if(!all(colnames(NPXbody)[2:93] == header[4,2:93])){
    message("ABORT ABORT ABORT", plateName)
  }
  
  tail <- plate %>% group_by(Assay) %>% 
    mutate(belowLOD = paste0(mean(NPX < LOD, na.rm = T)*100, "%"), normalization = "BAMBOO") %>% 
    distinct(Assay, LOD, belowLOD, normalization, .keep_all = T) %>% 
    select(Assay, LOD, belowLOD, normalization, AssayFlag) %>% 
    filter(LOD == max(LOD, na.rm = T)) %>% mutate(LOD = as.character(LOD)) %>% 
    arrange(match(Assay, header[4,2:93])) %>% t()
  
  if(!all(tail[1,] == header[4,2:93])){
    message("ABORT ABORT ABORT", plateName)
  }
  
  
  
  tail <- cbind(c(NA, "LOD", "Missing Data freq.", "Normalization", "below limit of detection adjustment"), tail)
  
  emptyMat <- matrix(nrow = 5, ncol = 4)
  
  tail <- cbind(tail,emptyMat)
  
  tail[1,] <- tail["AssayFlag",]
  tail["AssayFlag",] <- NA
  
  row.names(tail) <- NULL
  
  
  if(all(is.na(plate$`QC Deviation Det Ctrl`))){
    header <- header[,-c(96,97)]
    NPXbody <- NPXbody[,-c(96,97)]
    tail <- tail[,-c(96,97)]
  }
  
  top <- rbind(header, NPXbody)
  
  exportPlate <- rbind(top, tail)
  colnames(exportPlate) <- NULL
  
  if(!dir.exists(directory)){
    dir.create(directory)
  }
  
  file <- paste0(directory, plateName, "_", normalizationMethod, ".xlsx")
  
  openxlsx::write.xlsx(data.frame(exportPlate), file = file, colNames = F)
  read_NPX(file)
}

writeOlinkXLSX <- function(plateList, normMethod = "BAMBOO", directory){
  for(plateName in names(plateList)){
    plate <- plateList[[plateName]]
    headerAssayInfo <- plate %>% select(c(Panel, Assay, UniProt, OlinkID)) %>% unique() 
    header <- matrix(data = NA, nrow = 2, ncol = nrow(headerAssayInfo) + 5)
    header[1,2] <- "BAMBOO_export_v1"
    header[2,1] <- "NPX data"
    headerAssayInfo1 <- bind_rows(c("Panel" = "Panel", "Assay" = "Assay", "UniProt" = "Uniprot ID", "OlinkID" = "OlinkID"), headerAssayInfo) #%>% mutate(emptyRow = NA) #%>% t()
    headerAssayInfo2 <- data.frame(Panel = unique(plate$Panel),
                                   Assay = c("PlateID", "QC Warning", "QC Deviation from median", "QC Deviation from median"),
                                   UniProt = c(NA, NA, "Inc Ctrl 2", "Det Ctrl"),
                                   OlinkID = NA)
    headerAssayInfo <- bind_rows(headerAssayInfo1, headerAssayInfo2) %>% mutate(emptyRow = NA) %>% t()
    rownames(headerAssayInfo) <- NULL
    header <- rbind(header, headerAssayInfo)
    
    NPXbody <- plate %>% 
      mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% 
      select(c(SampleID_plate, NPX, Assay, PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`)) %>% 
      pivot_wider(names_from = c(Assay), values_from = NPX) %>% 
      relocate(c(PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`), .after = `CSF-1`) %>% 
      as.matrix()
    
    NPX <- plate %>% 
      mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% select(c(NPX, Assay, SampleID_plate)) %>% pivot_wider(names_from = Assay, values_from = NPX) %>% select(-SampleID_plate)
    
    tail <- plate %>%  dplyr::slice(c(seq(1,8096, 88))) %>% pull(LOD) 
    
    belowLOD <- apply(t(as.matrix(NPX)),1,function(x){paste0( round(( sum( x < tail ) / length( x ) ) * 100, digits = 2), "%")})
    
    tail <- cbind(matrix(c("LOD", tail, "Missing Data freq.", belowLOD, "Normalization", rep(normMethod, length(belowLOD))), nrow = 3, byrow = T), matrix(nrow = 3, ncol = 6))
    
    tail <- rbind(matrix(nrow = 1, ncol = ncol(tail)), tail)
    
    exportPlate <- rbind(header,NPXbody,tail)
    
    file <- paste0(directory, plateName, "_", normMethod, ".xlsx")
    
    openxlsx::write.xlsx(exportPlate, file = file, col.names = F)
  }
  
  # plates <- Reduce(bind_rows, plateList)
  # 
  # headerAssayInfo <- plates %>% select(c(Panel, Assay, UniProt, OlinkID)) %>% unique() 
  # header <- matrix(data = NA, nrow = 2, ncol = nrow(headerAssayInfo) + 5)
  # header[1,2] <- "Hidde_OlinkManager_0.1"
  # header[2,1] <- "NPX data"
  # headerAssayInfo1 <- bind_rows(c("Panel" = "Panel", "Assay" = "Assay", "UniProt" = "Uniprot ID", "OlinkID" = "OlinkID"), headerAssayInfo) #%>% mutate(emptyRow = NA) #%>% t()
  # headerAssayInfo2 <- data.frame(Panel = unique(plates$Panel),
  #                                Assay = c("Plate ID", "QC Warning", "QC Deviation from median", "QC Deviation from median"),
  #                                UniProt = c(NA, NA, "Inc Ctrl 2", "Det Ctrl"),
  #                                OlinkID = NA)
  # headerAssayInfo <- bind_rows(headerAssayInfo1, headerAssayInfo2) %>% mutate(emptyRow = NA) %>% t()
  # rownames(headerAssayInfo) <- NULL
  # header <- rbind(header, headerAssayInfo)
  # 
  # NPXbody <- plates %>% 
  #   mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% 
  #   select(c(SampleID_plate, NPX, Assay, PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`)) %>% 
  #   pivot_wider(names_from = c(Assay), values_from = NPX) %>% 
  #   relocate(c(PlateID, QC_Warning, `QC Deviation Inc Ctrl`, `QC Deviation Det Ctrl`), .after = `CSF-1`) %>% 
  #   as.matrix()
  # 
  # NPX <- plates %>% 
  #   mutate(SampleID_plate = paste0(SampleID, "_", plate)) %>% select(c(NPX, Assay, SampleID_plate)) %>% pivot_wider(names_from = Assay, values_from = NPX) %>% select(-SampleID_plate)
  # 
  # tail <- rep(NA, length(unique(plates$Assay)))
  # 
  # belowLOD <- rep(NA, length(unique(plates$Assay)))
  # 
  # tail <- cbind(matrix(c("LOD", tail, "Missing Data freq.", belowLOD, "Normalization", rep(normMethod, length(belowLOD))), nrow = 3, byrow = T), matrix(nrow = 3, ncol = 4))
  # 
  # tail <- rbind(matrix(nrow = 1, ncol = ncol(tail)), tail)
  # 
  # exportPlate <- rbind(header,NPXbody,tail)
  # 
  # file <- paste0(directory, plateName, "_AllPlates_", normMethod, ".xlsx")
  # 
  # openxlsx::write.xlsx(exportPlate, file = file, col.names = F)
}
