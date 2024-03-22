suppressWarnings(suppressMessages(require(tidyverse))) #
suppressWarnings(suppressMessages(require(openxlsx))) #
suppressWarnings(suppressMessages(require(readr))) #
suppressWarnings(suppressMessages(require(ggplot2))) #
suppressWarnings(suppressMessages(require(robustbase)))
suppressWarnings(suppressMessages(require(gridExtra))) #

recode <- dplyr::recode
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
group_by <- dplyr::group_by
case_when <- dplyr::case_when

# This scripts has the functions used to normalize Olink NPX data using 
# bridging controls. The functions are:
# - loadNPXfiles: loads all NPX files in a folder and returns a list of dataframes
# - readNPX: reads a single NPX file and returns a long format dataframe
# - BAMBOO_normalization: normalizes the data between two plates using bridging controls
# - BAMBOO_plot: plots the data before and after normalization
# - write_NPX_files: writes the normalized data in wide format to a folder

message("                                                                                
                           %%%%%%%%%%      .%%%%%%%%%&                          
                           #%%%%%%%%%#     %%%%%%%%%%,                          
                *%%%%.      %%%%%%%%%%/   %%%%%%%%%%%      ,%%%%.               
             ,%%%%%%%%%*    %%%%% .%%%%, *%%%%/ %%%%(    (%%%%%%%%%,            
            *%%%%%%%%%%%#(  *%%%%   %%%%.%%%%.  (%%%.  %%#%%%%%%%%%%,           
              #%%%%# *%%%%%%.%%%#    %%%%%%/    #%%%.%%%%%%,(%%%%%(             
                %%%%%/   #%%%%%%.     %%%#      %%%%%%%%*  .%%%%&               
     ,%%%%#,     ,%%%%/    .%%%%       (,       %%%#      ,%%%%.     *%#%%%.    
    (%%%%%%%%%%%%* (%%%.      .%                ,        #%%%/./%%%%%%%%%%%%/   
   #%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%/%%%%%%%(  
   *%%%%%%%%&,     ./#%%%                              (.           #%%%%%%%#.  
       .%%%%%%%%             %%%%%%,         *%%%%%%             (%%%%%%#       
            (%%%%%           %%%%%%,         *%%%%%%         *%%%%#&/           
 ,*/(#%%%&%%%%%%%%%%/        %%%%%%,         *%%%%%%       %%%%%%%%%%%%%&%%#(/**
*%%%%%%%%%%%%%%#/*.          %%%%%%,         *%%%%%%             (%%%%%%%%%%%%%%
/%%%%%%%,                    %%%%%%,         *%%%%%%                    /%%%%%%%
,%%%%%%%%%%%%%%(.            %%%%%%,         *%%%%%%         .*(#%%%%%%%%%%%%%%%
 ..,*/(#%&%%%%%%%%%%%#       %%%%%%*         /%%%%%#        *%%%%%%%%%%%%#(/*,..
            #%%%%%%.         *%%%%%%         %%%%%%,          .%%%%%#           
       ,%%%%%%%,              /%%%%%%%,   ,%%%%%%%/              #%%%%%%%,      
   /%%%%%%%%/          ./%      %#%%%%%%%%%%%%%%#       #%%%#.     .&%%%%%%%%/  
   #%%%%%%%/%%%%%%%%%%%%%          ,(%%%%%%#(.           %%%%%%%%%%%%%%%%%%%%(  
    (%%%%%%%%%%%#, #%%%(        *                %,      .%%%# ,#%%%%%%%%%%%/   
     ,%#%%(      ,%%%%      ,%%#%       /(       %%%%,    (%%%%,      #%%%%.    
               .%%%%%   /%%%%%%%#      %%%%     ,%%%%%%#   (%%%%%               
              #%%%%%*/%%%%%#.%%%(    #%%%%%#    %%%% %%%%%%(.#%%%%#             
            *%%%%%%%%%%%#(  ,%%%/  /%%%%.%%%%  .%%%%,  (%%%%%%%%%%%%*           
             ,%%%%%%%%%,    #%%%# #%%%%, ,%%%%. %%%%#    ,%%%%%%%%%.            
                .&%%%       %%%%%%%%%%/   #%%%%%%%%%%      .%%%%.               
                           (%%%%%%%%%%     %%%%%%%%%%*                          
                           %%%%%%%%%%      .%%%%%%%%%&                          

        
        
                        Importing the BAMBOO functions 
                          h.m.smits-7@umcutrecht.nl")

loadNPXfiles <- function(path){
  plateFiles <- list.files(path, full.names = T)
  plates <- Reduce(bind_rows, lapply(plateFiles, readNPX))
  plateList <- split(plates,f = plates$PlateID)
  return(plateList)
}
# file <- "/Users/Hidde/Documents/GitHub/BAMBOO/data//normalizedSubjectPlate.xlsx"  
readNPX <- function(file){
  suppressMessages(p <- read.xlsx(file, colNames = F))
  if(p[2,1] == "NPX data"){
    ########################################################################
    # This is classic style Olink Data, we transform it to long format data#
    ########################################################################
    header <- p[1:6,]
    EndBodyRow <- which(p[,1] == "LOD")
    NPXdata <- p[7:EndBodyRow-1,1:93]
    SampleData <- p[7:EndBodyRow-1,94:97]
    tail <- p[(EndBodyRow ):(EndBodyRow + 2),]
    
    Assays <- unname(unlist(header[4,2:93]))
    LOD_df <- data.frame(Assay = Assays, LOD = unname(unlist(tail[1,2:93])))
    plateID <- data.frame(SampleID = NPXdata[,1], 
                          PlateID = SampleData[,1],
                          stringsAsFactors = F
                          )
    
    colnames(NPXdata) <- c("SampleID", Assays)
    
    suppressMessages(
      longData <- NPXdata %>% 
      pivot_longer(cols = all_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
      mutate(NPX = as.numeric(NPX)) %>% 
      left_join(plateID) %>% 
      left_join(LOD_df) %>% 
      mutate(LOD = as.numeric(LOD)))
      
  }else{
    colnames(p) <- p[1,]
    p <- p[-1,]
    Assays <- colnames(p)[3:94]
    LOD <- data.frame(LOD = unlist(p[which(p[,1] == "LOD"),3:94] ), Assay = Assays)
    
    NPXdata <- p[1:82,] 
    suppressMessages(longData <- pivot_longer(NPXdata, cols = all_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
      left_join(LOD)  %>% 
      mutate(NPX = as.numeric(NPX)) %>% 
      mutate(LOD = as.numeric(LOD)))
    
  }
  return(longData)
}

BAMBOO_normalization <- function(plateReference, plateSubject, BCs, LODthreshold = 6){
  # The original BAMBOO function that normalizes the data between two plates using bridging controls
  # plateReference: the references plate as NPX long format file
  # plateSubject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # LODthreshold: Number of values that have to be above LOD to be not flagged
  
  mutate <- dplyr::mutate
  if(is.null(BCs)){
    BCs <- intersect(plateReference$SampleID, plateSubject$SampleID)
  }
  
  BCs <- removeOutliers(plateReference, plateSubject, BCs, quantileThreshold = 0.95) #based on interquantile outlier detection
  
  flaggedAssays <- flagAssay(plateReference = plateReference, plateSubject = plateSubject, BCs = BCs, BCsAboveLOD = LODthreshold, correlationThreshold = 0)
  
  plateReferenceBCs <- plateReference %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays) %>% 
    group_by(Assay) %>% 
    filter(SampleID%in%BCs)  %>% 
    mutate(NPX = case_when(NPX >= LOD ~ NPX, # check if order is correct
                           AssayFlag  ~ NPX,
                           NPX < LOD ~ NA))  
  

  
  plateSubjectBCs <- plateSubject %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays) %>% 
    group_by(Assay) %>% 
    filter(SampleID%in%BCs) %>% 
    mutate(NPX = case_when(NPX >= LOD ~ NPX, # check if order is correct
                           AssayFlag  ~ NPX,
                           NPX < LOD ~ NA)) 
    

  
  mergedPlates <- inner_join(plateReferenceBCs, plateSubjectBCs, by = c("SampleID", "Assay"), suffix = c(".ref", ".sub")) 
  
  if(!all(mergedPlates$SampleID%in%BCs)){
    
    message("sample ", setdiff(mergedPlates$SampleID, BCs), "not present in both plates")
    
  }
  
  modelCoeff <- lmrob(NPX.ref ~ NPX.sub, data = mergedPlates, na.action = "na.omit")$coefficients
  
  Adj_factors <- mergedPlates %>% 
    group_by(Assay) %>% 
    mutate(Adj_factor = median(NPX.ref - (NPX.sub*modelCoeff[["NPX.sub"]] + modelCoeff[["(Intercept)"]]), na.rm = T)) %>% 
    dplyr::select(c("Assay", "Adj_factor")) %>% 
    unique() %>% ungroup()
  
  plateSubject <- plateSubject %>% 
    left_join(Adj_factors, by = "Assay") %>% 
    dplyr::mutate(NPX = NPX*modelCoeff[["NPX.sub"]] + modelCoeff[["(Intercept)"]] + Adj_factor,
           LOD = LOD*modelCoeff[["NPX.sub"]] + modelCoeff[["(Intercept)"]] + Adj_factor) %>% 
    mutate(AssayFlag = Assay%in%flaggedAssays)
  
  return(plateSubject)
}



flagAssay <- function(plateReference, plateSubject, BCs, plateBelowLOD = 1, BCsAboveLOD = 6, correlationThreshold = 0.01){

  belowLODBCs.ref <- plateReference %>% filter(SampleID%in%BCs) %>% mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% group_by(Assay) %>% filter((sum(NPX > LOD) < BCsAboveLOD)) %>% pull(Assay) %>% unique()
  belowLODBCs.sub <- plateSubject %>% filter(SampleID%in%BCs) %>% mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% group_by(Assay) %>% filter((sum(NPX > LOD) < BCsAboveLOD)) %>% pull(Assay) %>% unique()
  
  # message("The following assays were removed because they were below detection limit in ", plateBelowLOD*100, "% of the samples ", paste(unique(c(belowLODplate.ref, belowLODplate.sub)), collapse = ", ") )
  message("The following assays were flagged because less than ", BCsAboveLOD, " samples were above detection limit ",paste(unique(c(belowLODBCs.ref, belowLODBCs.sub)), collapse = ", "))

  flaggedAssays <- unique(c(belowLODBCs.ref, belowLODBCs.sub))
  
  return(flaggedAssays)
}


removeOutliers <- function(plateReference, plateSubject, BCs, quantileThreshold = 0.95){ #inter Quantile outlier detection -> can be improved
  
  ss <- sumOfSquares(plateReference, plateSubject, samplesOfInterest = BCs, na.rm = T)
  
  outlierSamples <- ss %>% filter(SS > quantile(ss$SS, quantileThreshold)) %>% pull(SampleID)
  
  message(paste0("identiefied the following sample as outlier: ", outlierSamples, "\n"))
  
  return(BCs[!BCs%in%outlierSamples])
  
}

sumOfSquares <- function(plateRef, plateSubject, samplesOfInterest = NULL, na.rm = F){ 
  if(is.null(samplesOfInterest)){
    samplesOfInterest <- intersect(plateRef$SampleID, plateSubject$SampleID)
  }
  ss <- inner_join(plateRef %>% filter(SampleID %in% samplesOfInterest), plateSubject %>% filter(SampleID %in% samplesOfInterest), by = c("SampleID", "Assay")) %>% 
    filter(SampleID%in%samplesOfInterest) %>% group_by(SampleID) %>% 
    summarise(SS = sum((NPX.x - NPX.y)**2, na.rm = na.rm))
  return(ss)
}


writeNPX <- function(plate, path, filename){
  # This functions save the long format NPX data in a excel file
  
  body <- plate %>% select(-Adj_factor, -LOD, -AssayFlag, -plate) %>% pivot_wider(names_from = "Assay", values_from = "NPX")
  plateID <- plate %>% pull(plate) %>% unique()
  tail <- plate %>% 
    select(Adj_factor, LOD, AssayFlag) %>% 
    mutate(AssayFlag = AssayFlag * 1) %>% 
    distinct() %>% 
    t() %>% 
    as.data.frame() # this must be merged based on the assays, below the body
  blankDF <- data.frame(matrix(data = NA, nrow = nrow(tail), ncol = 2))%>% mutate_all(as.character) 
  tail <- cbind(blankDF, tail)[-4,] 
  colnames(tail) <- colnames(body)
  
  NPX_data <- bind_rows(body, tail) 
  NPX_data$plateID <- plateID
  write.xlsx(NPX_data, file = paste0(path, filename))
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
  

}

plotBeforeAndAfter <- function(referencePlate, subjectPlate, norm.SubjectPlate){
  
  p1 <- inner_join(referencePlate, subjectPlate, by = c("SampleID", "Assay")) %>% ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + geom_point() + theme_bw() + theme(legend.position = "none") + geom_abline(linetype = 2) + xlim(c(-2,16)) + ylim(c(-2,16))  + ggtitle("Before")
  p2 <- inner_join(referencePlate, norm.SubjectPlate, by = c("SampleID", "Assay")) %>% ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + geom_point() + theme_bw() + theme(legend.position = "none") + geom_abline(linetype = 2)+ xlim(c(-2,16))+ ylim(c(-2,16)) + ggtitle("After")
  
  gridExtra::grid.arrange(p1,p2, nrow = 1)
}

renameSamplesAndAddPlate <- function(data){
  ##### mostly bridging controls have multiple names on various plates, this renames them to the correct name
  ##### also renames BCs which where thawed mutliple times and looked to be good
  ##### Add a column with plate numbers 
  
  data <- data %>% mutate(
    SampleID = recode(
      SampleID,
      #plate 2 BC to plate 1
      "JDM#1" = "BC_JDM_p1", #
      "JDM#3" = "BC_JDM_p3", #
      "JDM#4" = "BC_JDM_p4", #
      "JDM#7 heparin" = "BC_JDM_p7", #
      "JDM#8 heparin" = "BC_JDM_p8", # 
      "JDM#9 heparin" = "BC_JDM_p9", #
      "HC_1 heparin" = "HC_Sohep_1", #
      "HC_7" = "HC-Sohep_7", #
      "IBD_7 heparin" = "Twin7_Sohep", #
      "HC_6 heparin" = "HC_Sohep_6", #
      "IBD_4 heparin" = "Twin4_Sohep", #
      "IBD_2" = "Twin2_Sohep", # 
      "JDM #5" = "BC_JDM_s5", #
      "JDM#7 serum" = "BC_JDM_s7", #
      "JDM#8 serum" = "BC_JDM_s8", #
      "JDM#9 serum" = "BC_JDM_s9", #
      "IDB_7" = "Twin7_serum" , #
      "HC_1 serum" = "HC_serum_1_zonder_gel", #
      "HC_6 serum" = "HC_serum_6_met_gel", #
      "IBD_1" = "Twin1_serum", #
      "IBD_6" = "Twin6_serum", #
      "HC_9" = "HC_serum_9_zonder_gel", #
      "HC_4" = "HC_serum_4_zonder_gel", #
      "IBD_4 serum" = "Twin4_serum", #
      "HC_2" = "HC_serum_2_met_gel", # 
      
      # different naming same samples because of replicates!
      "150903-55062 3-9-2015" = "onset_JDMA003_s",
      "JDM0287 22-1-2016" = "onset_ JDMU009_s",
      "P17-58232 18-5-2017" = "onset_JDMN004_s",
      "P17-59042 2-6-2017" = "onset_JDMR006_s",
      "P17-61987 20-7-2017" = "onset_ JDMU033_s",
      "P17-72116 21-12-2017" = "onset_JDMN005_s",
      "P18-66720 5-7-2018" = "active_JDMN008_s",
      "P18-76255 4-10-2018" = "active_JDMN006_s",
      "P19-51219 17-1-2019" = "remission_JDMN007_s",
      "P17-56441 19-4-2017" = "onset_JDMU027_2_p",
      
      #plate 3,4,5 to plate 1 BCs
      "JDM#4_p" = "BC_JDM_p4", #
      "JDM#7_p" = "BC_JDM_p7", #
      "JDM#9_p" = "BC_JDM_p9", #
      "HC_1_p" = "HC_Sohep_1", #
      "HC_7_p" = "HC-Sohep_7", #
      "IBD_2_p" = "Twin2_Sohep", #
      "JDM#5_s" = "BC_JDM_s5", #
      "JDM#8_s" = "BC_JDM_s8", #
      "HC_6_wg_s" = "HC_serum_6_met_gel", #
      "HC_9_wog_s" = "HC_serum_9_zonder_gel", #
      "IBD_4_s" = "Twin4_serum", #
      "HC_2_wg_s" = "HC_serum_2_met_gel", #
      "HC_6,_wg_s" = "HC_serum_6_met_gel",
      "IBD_4" = "Twin4_serum"
      
      
    )
  ) %>% group_by(PlateID) %>% 
    mutate(SampleID = case_when(
      SampleID == "IBD_4" ~ "Twin4_serum",
      SampleID == "JDM#7_1" & !"BC_JDM_p7"%in%SampleID ~ "BC_JDM_p7",
      # SampleID == "IBD_4_1"~ "Twin4_serum",
      # SampleID == "IBD_4_2"~ "Twin4_serum",
      SampleID == "IBD_2_1" & !"Twin2_Sohep"%in%SampleID ~ "Twin2_Sohep",
      SampleID == "JDM#4_2" & !"BC_JDM_p4"%in%SampleID ~ "BC_JDM_p4",
      SampleID == "JDM#7"~ "BC_JDM_p7",
      SampleID == "JDM#5"~ "BC_JDM_s5",
      # SampleID == "JDM#7_2" & !"BC_JDM_p7"%in%SampleID ~ "BC_JDM_p7",
      SampleID == "JDM#8_3" & !"BC_JDM_s8"%in%SampleID ~ "BC_JDM_s8",
      # SampleID == "JDM#9_2" & !"BC_JDM_p9"%in%SampleID ~ "BC_JDM_p9",
      SampleID == "IBD_4_1"& !"Twin4_serum"%in%SampleID ~ "Twin4_serum",
      # SampleID == "IBD_2_1"& !"Twin2_Sohep"%in%SampleID ~ "Twin2_Sohep",
      SampleID == "JDM#5_2" & !"BC_JDM_s5"%in%SampleID ~ "BC_JDM_s5",
      # SampleID == "JDM#7_2" & !"BC_JDM_s7"%in%SampleID ~ "BC_JDM_s7",
      # SampleID == "JDM#8_s_1/27" & !"BC_JDM_s8"%in%SampleID ~ "BC_JDM_s8",
      # SampleID == "JDM#9_2" & !"BC_JDM_s9"%in%SampleID ~ "BC_JDM_s9",
      # SampleID == "JDM#4,_2" & !"BC_JDM_p4"%in%SampleID ~ "BC_JDM_p4",
      # SampleID == "JDM#8_s_1/27" & !"BC_JDM_p8"%in%SampleID ~ "BC_JDM_p8",
      SampleID == "IBD_4_2" & !"Twin4_serum"%in%SampleID & ! "IBD_4_1" %in% SampleID  ~ "Twin4_serum", #& !"IBD_4_1"%in%SampleID & !"IBD_4"%in%SampleID
      SampleID == "IBD_4_4" & !"Twin4_serum"%in%SampleID ~ "Twin4_serum", #& !"IBD_4_1"%in%SampleID & !"IBD_4"%in%SampleID 
      SampleID == "IBD_2_2" & !"Twin2_Sohep"%in%SampleID & !"IBD_2_1"%in%SampleID  ~ "Twin2_Sohep" , #& !"IBD_2_1"%in%SampleID
      SampleID == "HC_2, met gel" ~ "HC_serum_2_met_gel",
      SampleID == "HC_6, met gel" ~ "HC_serum_6_met_gel",
      SampleID == "HC_9 zonder gel" ~ "HC_serum_9_zonder_gel",
      SampleID == "HC_1" ~ "HC_Sohep_1",
      SampleID == "JDM#4, JDMN006" ~ "BC_JDM_p4",
      SampleID == "JDM#9, JDML001" ~ "BC_JDM_p9",
      SampleID == "JDM#7, JDMN010" ~ "BC_JDM_p7",
      SampleID == "JDM#8, JDMR008" ~ "BC_JDM_s8",
      SampleID == "JDM#5, JDMN007" ~ "BC_JDM_s5",
      SampleID == "JDM#8_s_1/27" ~ "BC_JDM_s8",
      SampleID == "JDM#4, aliquot 2e batch" & !"JDM#4, JDMN006" %in% SampleID & ! "BC_JDM_p4" %in% SampleID~ "BC_JDM_p4",
      SampleID == "JDM#9" ~ "BC_JDM_p9",
      SampleID == "JDM#4,_2" & !"BC_JDM_p4"%in%SampleID  ~ "BC_JDM_p4",
      # SampleID == "IBD_4_1" ~ "Twin4_serum",
      # SampleID == "JDM#8_3" ~ "BC_JDM_s8",
      # SampleID == "IBD_2_2" ~ "Twin2_Sohep",
      # SampleID == "JDM#4,_2" ~ "BC_JDM_p4",
      SampleID == "JDM#8_3" & !"BC_JDM_s8"%in%SampleID ~ "BC_JDM_s8",
      TRUE ~ SampleID
    ))%>% 
    mutate(plate = 
             case_when(
               PlateID == "19aug2021_Zo project_plt 15"~15,
               PlateID == "19may_SN&ED_pilot#4_plt6 I-O"  ~ 6,
               PlateID == "20may2021_SN&ED_pilot#34_plt7 I-O" ~ 7,
               PlateID == "20may2021_SN&ED_pliot#4_plt8 I-O"    ~ 8,
               PlateID == "23apr2021_SN&ED_pilot#3_plt3_I-O"    ~ 3,
               PlateID == "23apr2021_SN&ED_pilot#3_plt4_I-O"   ~ 4,
               PlateID == "23apr2021_SN&ED_pilot#3_plt5_I-O"  ~ 5,
               PlateID ==  "24aug2021_Zo project_plt 16"   ~ 16,
               PlateID == "24mrt2021_SNierkens&EDelemarre_pilot#2_I-O"  ~ 2,
               PlateID == "26feb2021_Stefan Nierkens_I-O_pilot"   ~ 1,
               PlateID == "7jun2021_Zo project_I-O_plt 10" ~ 10,
               PlateID ==  "7jun2021_Zo project_I-O_plt 11"  ~ 11,
               PlateID ==  "7jun2021_Zo project_I-O_plt 9" ~ 9,
               PlateID ==  "mrt2022_Zo_I-O_plt 25"   ~ 25,
               PlateID ==  "mrt2022_Zo_I-O_plt 28"     ~ 28,
               PlateID == "mrt2022_Zo_I-O_plt 34" ~ 34,
               PlateID == "mrt2022_Zo_I-O_plt 36"  ~ 36,
               PlateID == "mrt2022_Zo-I-O_plt 29"  ~ 29,
               PlateID == "mrt2022-Zo_I-O_plt 35"  ~ 35,
               PlateID == "Zo project_I-O_plt 26"  ~ 26,
               PlateID == "Zo project_I-O-plt 27" ~ 27,
               PlateID == "Zo project_jun2021_plt 12" ~ 12,
               PlateID == "Zo project_jun2021_Plt 13" ~ 13,
               PlateID == "Zo project_jun2021_plt 14" ~ 14,
               PlateID == "Zo_plt21_I-O"  ~ 21,
               PlateID == "Zo_plt22_I-O" ~ 22,
               PlateID ==  "ZO_plt32_I-O"  ~ 32,
               PlateID == "ZO_plt33_I-O" ~ 33,
               PlateID == "Zo_project_plt 23" ~ 23,
               PlateID == "Zo-plt 17" ~ 17,
               PlateID ==  "Zo-plt 18" ~ 18,
               PlateID ==  "Zo-plt 19" ~ 19,
               PlateID == "Zo-plt 20"  ~ 20,
               PlateID == "Zo-project_plt 24" ~ 24,
               PlateID == "26apr2022_Zo_I-O_plt 37" ~ 37,
               PlateID == "26apr2022_Zo_I-O_plt 38" ~ 38,
               PlateID == "26apr2022_Zo_O-D_plt 37" ~ 37,
               PlateID == "26apr2022_Zo_O-D_plt 38" ~ 38,
               PlateID == "ZO_plt31_OD" ~ 31,
               PlateID == "Zo_plt32_OD" ~ 32, 
               PlateID ==  "Zo_plt33_OD" ~33,
               PlateID == "mrt2022_Zo_O-D_plt 34" ~ 34,
               PlateID == "mrt2022-Zo_O-D_plt 35" ~ 35,
               PlateID == "Zo_O-D_plt 36" ~ 36,
               PlateID == "13jul2022_Zo_I-O_plt 30" ~ 30,
               PlateID == "7jul2022_Zo_I-O_plt 32" ~ 32,
               PlateID == "7jul2022_Zo_I-O_plt 39" ~ 39,
               PlateID == "14jul2022_Zo_I-O_plt 31" ~ 31,
               PlateID == "7jul2022_Zo_O-D_plt 39" ~ 39,
               PlateID == "Zo project_I-O_plt 40" ~ 40,
               PlateID == "Zo project_O-D_plt 40" ~ 40,
               T ~ NA_real_
             )) %>% ungroup()
  return(data)
}



