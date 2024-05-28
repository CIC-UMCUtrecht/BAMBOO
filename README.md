# BAMBOO Batch Correction

This repository contains the functions for performing BAMBOO batch correction as described in the accompanying paper "Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method."

The BAMBOO batch correction method is designed to correct batch effects in Olink Explore data. It is based on a robust linear regression model that corrects for the three types of batch effects that were identified (plate, sample and protein) whilst presevering biological variability.

## Usage
To use BAMBOO batch correction easily done in the R environment. To use BAMBOO batch correction in your own project simply make a vector with the names of the bridging controls in your data and run the BAMBOO function. The function will return a corrected data frame with the batch effects removed.

Example Script: Check out the example script provided in this repository to see how to apply BAMBOO batch correction to your Olink explore data.

## output
- LOD
- Files
- plots

## Citation
If you find BAMBOO batch correction useful in your research, please consider citing our paper:
"Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method."
- DOI
