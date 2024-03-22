#####################################################
################ Example script #####################
#####################################################

source("/Users/Hidde/Documents/GitHub/BAMBOO/BAMBOO_functions.R")
path <- "/Users/Hidde/Documents/GitHub/BAMBOO/data/"
plateList <- loadNPXfiles("/Users/Hidde/Documents/GitHub/BAMBOO/data/")

referencePlate <- plateList[[3]] %>% renameSamplesAndAddPlate()
subjectPlate <- plateList[[1]] %>% renameSamplesAndAddPlate()

BCs <- intersect(referencePlate$SampleID, subjectPlate$SampleID)

norm.SubjectPlate <- BAMBOO_normalization(referencePlate, subjectPlate, BCs = BCs, LODthreshold = 6)

plotBeforeAndAfter(referencePlate, subjectPlate, norm.SubjectPlate)

writeNPX(norm.SubjectPlate, path = path, filename = "normalizedSubjectPlate.xlsx")
