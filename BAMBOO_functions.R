# Loading necessary packages
suppressWarnings(suppressMessages(require(openxlsx))) # Only package that is really needed
suppressWarnings(suppressMessages(require(robustbase))) # Only package that is really needed
suppressWarnings(suppressMessages(require(gridExtra))) # save plots separate

# Define dplyr functions for convenience
recode <- dplyr::recode
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
group_by <- dplyr::group_by
case_when <- dplyr::case_when

# This script has the functions used to normalize Olink NPX data using bridging controls. 
# The functions are:
# - load_NPX_files: loads all NPX files in a folder and returns a list of dataframes
# - read_NPX: reads a single NPX file and returns a long format dataframe
# - BAMBOO_normalization: normalizes the data between two plates using bridging controls
# - BAMBOO_plot: plots the data before and after normalization
# - write_NPX: writes the normalized data in wide format to a folder

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

load_NPX_files <- function(path) {
  # This function loads all NPX files in a folder and returns a list of dataframes
  # param path: path to the folder with the NPX files
  plate_files <- list.files(path, full.names = TRUE)
  plates <- Reduce(bind_rows, lapply(plate_files, read_NPX))
  plate_list <- split(plates, f = plates$PlateID)
  return(plate_list)
}

read_NPX <- function(file) {
  # This function reads a single NPX file and returns a long format dataframe
  # param file: path to the NPX file
  
  suppressMessages(p <- read.xlsx(file, colNames = FALSE))
  
  if (p[2, 1] == "NPX data") {
    # Classic style Olink Data, transform it to long format data
    header <- p[1:6,]
    EndBodyRow <- which(p[, 1] == "LOD")
    NPXdata <- p[7:(EndBodyRow - 1), 1:93]
    SampleData <- p[7:(EndBodyRow - 1), 94:97]
    tail <- p[EndBodyRow:(EndBodyRow + 2),]
    
    Assays <- unname(unlist(header[4, 2:93]))
    LOD_df <- data.frame(Assay = Assays, LOD = unname(unlist(tail[1, 2:93])))
    plateID <- data.frame(SampleID = NPXdata[, 1], PlateID = SampleData[, 1], stringsAsFactors = FALSE)
    
    colnames(NPXdata) <- c("SampleID", Assays)
    
    suppressMessages(
      longData <- NPXdata %>% 
        pivot_longer(cols = any_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
        mutate(NPX = as.numeric(NPX)) %>% 
        left_join(plateID) %>% 
        left_join(LOD_df) %>% 
        mutate(LOD = as.numeric(LOD))
    )
      
  } else {
    colnames(p) <- p[1, ]
    p <- p[-1, ]
    Assays <- colnames(p)[3:94]
    LOD <- data.frame(LOD = unlist(p[which(p[, 1] == "LOD"), 3:94]), Assay = Assays)
    
    NPXdata <- p[1:82, ] 
    suppressMessages(longData <- pivot_longer(NPXdata, cols = any_of(Assays), values_to = "NPX", names_to = "Assay") %>% 
      left_join(LOD)  %>% 
      mutate(NPX = as.numeric(NPX)) %>% 
      mutate(LOD = as.numeric(LOD)))
  }
  
  return(longData)
}

BAMBOO_normalization <- function(plate_reference, plate_subject, BCs, LOD_threshold = 6) {
  # The original BAMBOO function that normalizes the data between two plates using bridging controls
  # plate_reference: the references plate as NPX long format file
  # plate_subject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # LOD_threshold: Number of values that have to be above LOD to be not flagged
  
  if (is.null(BCs)) {
    BCs <- intersect(plate_reference$SampleID, plate_subject$SampleID)
  }
  
  BCs <- remove_outliers(plate_reference, plate_subject, BCs, quantile_threshold = 0.95) # based on interquantile outlier detection
  
  flagged_assays <- flag_assay(plate_reference, plate_subject, BCs, LOD_threshold)
  
  plate_reference_BCs <- plate_reference %>% 
    mutate(AssayFlag = Assay %in% flagged_assays) %>% 
    group_by(Assay) %>% 
    filter(SampleID %in% BCs)  %>% 
    mutate(NPX = case_when(NPX >= LOD ~ NPX, AssayFlag ~ NPX, NPX < LOD ~ NPX))
  
  plate_subject_BCs <- plate_subject %>% 
    mutate(AssayFlag = Assay %in% flagged_assays) %>% 
    group_by(Assay) %>% 
    filter(SampleID %in% BCs) %>% 
    mutate(NPX = case_when(NPX >= LOD ~ NPX, AssayFlag ~ NPX, NPX < LOD ~ NPX))
  
  merged_plates <- inner_join(plate_reference_BCs, plate_subject_BCs, by = c("SampleID", "Assay"), suffix = c(".ref", ".sub")) 
  
  if (!all(merged_plates$SampleID %in% BCs)) {
    message("Sample ", setdiff(merged_plates$SampleID, BCs), " not present in both plates")
  }
  
  model_coeff <- lmrob(NPX.ref ~ NPX.sub, data = merged_plates, na.action = "na.omit")$coefficients
  
  Adj_factors <- merged_plates %>% 
    group_by(Assay) %>% 
    mutate(Adj_factor = median(NPX.ref - (NPX.sub * model_coeff[["NPX.sub"]] + model_coeff[["(Intercept)"]]), na.rm = TRUE)) %>% 
    dplyr::select(c("Assay", "Adj_factor")) %>% 
    unique() %>% 
    ungroup()
  
  plate_subject <- plate_subject %>% 
    left_join(Adj_factors, by = "Assay") %>% 
    dplyr::mutate(NPX = NPX * model_coeff[["NPX.sub"]] + model_coeff[["(Intercept)"]] + Adj_factor,
                  LOD = LOD * model_coeff[["NPX.sub"]] + model_coeff[["(Intercept)"]] + Adj_factor) %>% 
    mutate(AssayFlag = Assay %in% flagged_assays)
  
  return(plate_subject)
}

flag_assay <- function(plate_reference, plate_subject, BCs, BCs_above_LOD = 6) {
  # Function that flags assays that are below the detection limit in a certain number of the samples
  # plate_reference: the references plate as NPX long format file
  # plate_subject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # BCs_above_LOD: Number of values that have to be above LOD to be not flagged
  
  below_LOD_BCs.ref <- plate_reference %>% filter(SampleID %in% BCs) %>% 
    mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% 
    group_by(Assay) %>% 
    filter((sum(NPX > LOD) < BCs_above_LOD)) %>% pull(Assay) %>% unique()
  
  below_LOD_BCs.sub <- plate_subject %>% filter(SampleID %in% BCs) %>% 
    mutate_at(vars(NPX), ~replace(., is.na(.), 0)) %>% 
    group_by(Assay) %>% 
    filter((sum(NPX > LOD) < BCs_above_LOD)) %>% pull(Assay) %>% unique()
  
  message("The following assays were flagged because less than ", BCs_above_LOD, " samples were above detection limit ", paste(unique(c(below_LOD_BCs.ref, below_LOD_BCs.sub)), collapse = ", "))
  
  flagged_assays <- unique(c(below_LOD_BCs.ref, below_LOD_BCs.sub))
  
  return(flagged_assays)
}

remove_outliers <- function(plate_reference, plate_subject, BCs, quantile_threshold = 0.95) { 
  # Function that removes outliers based on the sum of squares of the difference between the reference and the subject plate
  # plate_reference: the references plate as NPX long format file
  # plate_subject: the plate that has to be normalized to the reference plate, in NPX long format
  # BCs: vector with BC names (as seen in the SampleID column)
  # quantile_threshold: the threshold for the sum of squares to be considered an outlier
  
  ss <- sum_of_squares(plate_reference, plate_subject, samples_of_interest = BCs, na.rm = TRUE)
  
  outlier_samples <- ss %>% filter(SS > quantile(ss$SS, quantile_threshold)) %>% pull(SampleID)
  
  message(paste0("Identified the following sample as outlier: ", outlier_samples, "\n"))
  
  return(BCs[!BCs %in% outlier_samples])
}

sum_of_squares <- function(plate_reference, plate_subject, samples_of_interest = NULL, na.rm = FALSE) { 
  # Function that calculates the sum of squares of the difference between the reference and the subject plate BCs
  # plate_reference: the references plate as NPX long format file
  # plate_subject: the plate that has to be normalized to the reference plate, in NPX long format
  # samples_of_interest: vector with BC names (as seen in the SampleID column)
  # na.rm: remove NA values
  
  if (is.null(samples_of_interest)) {
    samples_of_interest <- intersect(plate_reference$SampleID, plate_subject$SampleID)
  }
  ss <- inner_join(plate_reference %>% filter(SampleID %in% samples_of_interest), plate_subject %>% filter(SampleID %in% samples_of_interest), by = c("SampleID", "Assay")) %>% 
    filter(SampleID %in% samples_of_interest) %>% 
    group_by(SampleID) %>% 
    summarise(SS = sum((NPX.x - NPX.y)^2, na.rm = na.rm))
  return(ss)
}

write_NPX <- function(plate.DF, path, filename) {
  # This function saves the long format NPX data in an Excel file
  # plate.DF: the NPX data in long format
  # path: the path where the file has to be saved
  # filename: the name of the file
  
  body <- plate.DF %>% select(-any_of(c("Adj_factor", "LOD", "AssayFlag", "PlateID"))) %>% pivot_wider(names_from = "Assay", values_from = "NPX")
  plateID <- plate.DF %>% pull(PlateID) %>% unique()
  tail <- plate.DF %>% 
    select(Adj_factor, LOD, AssayFlag) %>% 
    mutate(AssayFlag = as.numeric(AssayFlag)) %>% 
    distinct() %>% 
    t() %>% 
    as.data.frame()
  blankDF <- data.frame(matrix(data = NA, nrow = nrow(tail), ncol = 2)) %>% mutate_all(as.character) 
  
  NPX_data <- bind_rows(body)
  NPX_data$PlateID <- plateID
  write.xlsx(NPX_data, file = paste0(path, filename))
}

plot_before_and_after <- function(reference_plate, subject_plate, norm_subject_plate) {
  # Function that plots the NPX values of the reference plate, the subject plate, and the normalized subject plate
  # reference_plate: the references plate as NPX long format file
  # subject_plate: the plate that has to be normalized to the reference plate, in NPX long format
  # norm_subject_plate: the normalized values of the subject plate
  
  p1 <- inner_join(reference_plate, subject_plate, by = c("SampleID", "Assay")) %>% 
    ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + 
    geom_point() + 
    theme_bw() + 
    theme(legend.position = "none") + 
    geom_abline(linetype = 2) + 
    xlim(c(-2, 16)) + 
    ylim(c(-2, 16)) + 
    ggtitle("Before")
  
  p2 <- inner_join(reference_plate, norm_subject_plate, by = c("SampleID", "Assay")) %>% 
    ggplot(aes(x = NPX.x, y = NPX.y, col = Assay, shape = NPX.x < LOD.x | NPX.y < LOD.y)) + 
    geom_point() + 
    theme_bw() + 
    theme(legend.position = "none") + 
    geom_abline(linetype = 2) + 
    xlim(c(-2, 16)) + 
    ylim(c(-2, 16)) + 
    ggtitle("After")
  
  gridExtra::grid.arrange(p1, p2, nrow = 1)
}
