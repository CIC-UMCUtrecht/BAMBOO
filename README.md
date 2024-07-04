# BAMBOO Batch Correction

This repository contains functions for performing BAMBOO batch correction as described in the paper "Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method." [https://doi.org/10.21203/rs.3.rs-4044125/v1]

The BAMBOO batch correction method is designed to correct batch effects in Olink Explore data. It utilizes a robust linear regression model to correct for the three identified types of batch effects (plate, sample, and protein) while preserving biological variability.

## Usage
BAMBOO batch correction can be easily performed in the R environment. To use BAMBOO batch correction in your project, create a vector with the names of the bridging controls in your data and run the BAMBOO function. The function will return a dataframe of the subject plate with the batch effects removed. Check out the example script provided in this repository to see how to apply BAMBOO batch correction to your Olink Explore data.

## Output
The output of the function is a data frame containing the batch-corrected data of the subject plate. Of note here is that the LOD chosen in this data frame is also corrected, according to the proteins correction factor. If you want to output multiple plates in one dataframe, we recommend using the highest LOD value. Lastly, a function is provided that plots the bridging control samples of the subject and the reference plate against each other before and after BAMBOO batch correction.

## Citation
If you find BAMBOO batch correction useful in your research, please consider citing our paper:
"Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method." [https://doi.org/10.21203/rs.3.rs-4044125/v1]
