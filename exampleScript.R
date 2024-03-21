#####################################################
################ Example script #####################
#####################################################

source("/Users/Hidde/Documents/GitHub/BAMBOO/BAMBOO_functions.R")
path <- "/Users/Hidde/Documents/GitHub/BAMBOO/data/"
plateList <- loadNPXfiles("/Users/Hidde/Documents/GitHub/BAMBOO/data/")

referencePlate <- plateList[[2]]
subjectPlate <- plateList[[3]]

BCs <- intersect(referencePlate$SampleID, subjectPlate$SampleID)

norm.SubjectPlate <- BAMBOO_normalization(referencePlate, subjectPlate, BCs = BCs, LODthreshold = 6)

