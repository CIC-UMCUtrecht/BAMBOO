# BAMBOO Batch Correction

This repository contains functions for performing BAMBOO batch correction as described in the paper "Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method." [https://doi.org/10.21203/rs.3.rs-4044125/v1]

The BAMBOO batch correction method is designed to correct batch effects in Olink Explore data. It utilizes a robust linear regression model to correct for the three identified types of batch effects (plate, sample, and protein) while preserving biological variability.

## Usage
BAMBOO batch correction can be easily performed in the R environment. To use BAMBOO batch correction in your project, create a vector with the names of the bridging controls in your data and run the BAMBOO function. The function will return a dataframe of the subject plate with the batch effects removed. Check out the example script provided in this repository to see how to apply BAMBOO batch correction to your Olink Explore data.

## Output
The output of the function is a data frame containing the batch-corrected data of the subject plate. Note that the LOD chosen in this data frame is corrected for the correction factor. If multiple plates are merged into one data frame, we recommend using the highest LOD value. The source file also includes a function that plots the bridging controls against each other before and after batch correction. If noticeable batch effects are present, this will be evident in the first plot as many points deviating from the first diagonal.

## Citation
If you find BAMBOO batch correction useful in your research, please consider citing our paper:
"Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method." [https://doi.org/10.21203/rs.3.rs-4044125/v1]
