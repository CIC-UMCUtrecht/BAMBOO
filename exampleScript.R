#####################################################
################ Example Script #####################
#####################################################

# Load necessary functions
source("./GitHub/BAMBOO/BAMBOO_functions.R")
data_path <- "./GitHub/BAMBOO/data/"

# Load NPX data files (only wide format is supported)
plate_list <- loadNPXfiles(data_path)

# Define reference and subject plates
reference_plate <- plate_list[[1]]
subject_plate <- plate_list[[2]]

# Define bridging controls
bridging_controls <- intersect(reference_plate$SampleID, subject_plate$SampleID)

# Perform BAMBOO normalization
normalized_subject_plate <- BAMBOO_normalization(
  reference_plate, 
  subject_plate, 
  BCs = bridging_controls, 
  LODthreshold = 6
)

# Plot bridging controls before and after BAMBOO normalization
plotBeforeAndAfter(reference_plate, subject_plate, normalized_subject_plate)

# Save the normalized data in wide format
writeNPX(norm.SubjectPlate, path =  "./GitHub/BAMBOO/normData/", filename = "normalizedSubjectPlate.xlsx")

