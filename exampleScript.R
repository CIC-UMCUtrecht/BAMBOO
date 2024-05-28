#####################################################
################ Example script #####################
#####################################################

source("./GitHub/BAMBOO/BAMBOO_functions.R")
path <- "./GitHub/BAMBOO/data/"

# Here you can enter the path to the NPX data, currently BAMBOO only support the wide format of the data.

plateList <- loadNPXfiles("./GitHub/BAMBOO/data/")

# Here you define both plates

referencePlate <- plateList[[1]]

# Here you define the plate you'd like to normalize to the reference plate

subjectPlate <- plateList[[2]]

# Here you can define the bridging controls

BCs <- intersect(referencePlate$SampleID, subjectPlate$SampleID) 

# BAMBOO

norm.SubjectPlate <- BAMBOO_normalization(referencePlate, subjectPlate, BCs = BCs, LODthreshold = 6)

# Plot the bridging controls on both plates against each other before and after BAMBOO.

plotBeforeAndAfter(referencePlate, subjectPlate, norm.SubjectPlate)

# Save the data in wide format

writeNPX(norm.SubjectPlate, path =  "./GitHub/BAMBOO/normData/", filename = "normalizedSubjectPlate.xlsx")
