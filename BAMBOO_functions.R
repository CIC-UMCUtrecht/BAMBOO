suppressWarnings(suppressMessages(require(tidyverse))) # rewrite without the tidyverse????
suppressWarnings(suppressMessages(require(openxlsx))) # Only package that is really needed
suppressWarnings(suppressMessages(require(robustbase))) # Only package that is really needed
suppressWarnings(suppressMessages(require(gridExtra))) # save plots separate????

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
                          cic@umcutrecht.nl")

loadNPXfiles <- function(path){
  # This function loads all NPX files in a folder and returns a list of dataframes
  # param path: path to the folder with the NPX files
  plateFiles <- list.files(path, full.names = T)
  plates <- Reduce(bind_rows, lapply(plateFiles, readNPX))
  plateList <- split(plates,f = plates$PlateID)
  return(plateList)
}

readNPX <- function(file){
  # To do: make the format for the NPX file more flexible
  # This function reads a single NPX file and returns a long format dataframe
  # param file: path to the NPX file this can be in Olink format or a data.frame with samples as rows, and columns as proteins
  # To do: import multiple version of the data
  
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
      pivot_longer(cols = any_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
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
    suppressMessages(longData <- pivot_longer(NPXdata, cols = any_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
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
  
  flaggedAssays <- flagAssay(plateReference = plateReference, plateSubject = plateSubject, BCs = BCs, BCsAboveLOD = LODthreshold)
  
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



flagAssay <- function(plateReference, plateSubject, BCs, BCsAboveLOD = 6){
  # Function that flags assays that are below the detection limit in a certain number of the samples
  # plateReference: the references plate as NPX long format file
  # plateSubject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # BCsAboveLOD: Number of values that have to be above LOD to be not flagged
  
  belowLODBCs.ref <- plateReference %>% filter(SampleID%in%BCs) %>% mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% group_by(Assay) %>% filter((sum(NPX > LOD) < BCsAboveLOD)) %>% pull(Assay) %>% unique()
  belowLODBCs.sub <- plateSubject %>% filter(SampleID%in%BCs) %>% mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% group_by(Assay) %>% filter((sum(NPX > LOD) < BCsAboveLOD)) %>% pull(Assay) %>% unique()
  
  # message("The following assays were removed because they were below detection limit in ", plateBelowLOD*100, "% of the samples ", paste(unique(c(belowLODplate.ref, belowLODplate.sub)), collapse = ", ") )
  message("The following assays were flagged because less than ", BCsAboveLOD, " samples were above detection limit ",paste(unique(c(belowLODBCs.ref, belowLODBCs.sub)), collapse = ", "))

  flaggedAssays <- unique(c(belowLODBCs.ref, belowLODBCs.sub))
  
  return(flaggedAssays)
}


removeOutliers <- function(plateReference, plateSubject, BCs, quantileThreshold = 0.95){ #inter Quantile outlier detection -> can be improved
  # Function that removes outliers based on the sum of squares of the difference between the reference and the subject plate
  # plateReference: the references plate as NPX long format file
  # plateSubject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # quantileThreshold: the threshold for the sum of squares to be considered an outlier
  
  ss <- sumOfSquares(plateReference, plateSubject, samplesOfInterest = BCs, na.rm = T)
  
  outlierSamples <- ss %>% filter(SS > quantile(ss$SS, quantileThreshold)) %>% pull(SampleID)
  
  message(paste0("identiefied the following sample as outlier: ", outlierSamples, "\n"))
  
  return(BCs[!BCs%in%outlierSamples])
  
}

sumOfSquares <- function(plateRef, plateSubject, samplesOfInterest = NULL, na.rm = F){ 
  # Function that calculates the sum of squares of the difference between the reference and the subject plate BCs
  # plateReference: the references plate as NPX long format file
  # plateSubject: the plate that has to be normalized to the reference plate, in NPX long format
  # samplesOfInterest: vector with BC names (as seen in the SampleID column)
  # na.rm: remove NA values
  
  if(is.null(samplesOfInterest)){
    samplesOfInterest <- intersect(plateRef$SampleID, plateSubject$SampleID)
  }
  ss <- inner_join(plateRef %>% filter(SampleID %in% samplesOfInterest), plateSubject %>% filter(SampleID %in% samplesOfInterest), by = c("SampleID", "Assay")) %>% 
    filter(SampleID%in%samplesOfInterest) %>% group_by(SampleID) %>% 
    summarise(SS = sum((NPX.x - NPX.y)**2, na.rm = na.rm))
  return(ss)
}


writeNPX <- function(plateDF, path, filename){
  # This functions save the long format NPX data in a excel file
  # plate: the NPX data in long format
  # path: the path where the file has to be saved
  # filename: the name of the file
  # To do: Make the output in the correct format so it can be used with the Olink Analyze package
  
  body <- plateDF %>% select(-any_of(c("Adj_factor", "LOD", "AssayFlag", "PlateID"))) %>% pivot_wider(names_from = "Assay", values_from = "NPX")
  plateID <- plateDF %>% pull(PlateID) %>% unique()
  tail <- plateDF %>% 
    select(Adj_factor, LOD, AssayFlag) %>% 
    mutate(AssayFlag = AssayFlag * 1) %>% 
    distinct() %>% 
    t() %>% 
    as.data.frame() # this must be merged based on the assays, below the body
  blankDF <- data.frame(matrix(data = NA, nrow = nrow(tail), ncol = 2))%>% mutate_all(as.character) 
  # tail <- cbind(blankDF, tail)[-4,] 
  # colnames(tail) <- colnames(body)
  
  NPX_data <- bind_rows(body) #tail 
  NPX_data$plateID <- plateID
  write.xlsx(NPX_data, file = paste0(path, filename))
}

plotBeforeAndAfter <- function(referencePlate, subjectPlate, norm.SubjectPlate){
  # Function that plots the NPX values of the reference plate, the subject plate and the normalized subject plate
  # referencePlate: the references plate as NPX long format file
  # subjectPlate: the plate that has to be normalized to the reference plate, in NPX long format
  # norm.SubjectPlate: the normalized values of the subject plate
  
  p1 <- inner_join(referencePlate, subjectPlate, by = c("SampleID", "Assay")) %>% ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + geom_point() + theme_bw() + theme(legend.position = "none") + geom_abline(linetype = 2) + xlim(c(-2,16)) + ylim(c(-2,16))  + ggtitle("Before")
  p2 <- inner_join(referencePlate, norm.SubjectPlate, by = c("SampleID", "Assay")) %>% ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + geom_point() + theme_bw() + theme(legend.position = "none") + geom_abline(linetype = 2)+ xlim(c(-2,16))+ ylim(c(-2,16)) + ggtitle("After")
  
  gridExtra::grid.arrange(p1,p2, nrow = 1)
}
