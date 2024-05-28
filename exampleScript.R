#####################################################
################ Example script #####################
#####################################################

source("./GitHub/BAMBOO/BAMBOO_functions.R")
path <- "./GitHub/BAMBOO/data/"
plateList <- loadNPXfiles("./GitHub/BAMBOO/data/")

referencePlate <- plateList[[1]]
subjectPlate <- plateList[[2]]

BCs <- intersect(referencePlate$SampleID, subjectPlate$SampleID)

norm.SubjectPlate <- BAMBOO_normalization(referencePlate, subjectPlate, BCs = BCs, LODthreshold = 6)

plotBeforeAndAfter(referencePlate, subjectPlate, norm.SubjectPlate)

writeNPX(norm.SubjectPlate, path =  "./GitHub/BAMBOO/normData/", filename = "normalizedSubjectPlate.xlsx")
