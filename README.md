# BAMBOO

BAMBOO Batch Correction

This repository contains the source code for performing BAMBOO batch correction as described in the accompanying paper "Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method.""

The BAMBOO batch correction method is designed to correct batch effects in Olink Explore data. It is based on a robust linear regression model that corrects for the three types of batch effects that we identified (plate, sample and protein) whilst presevering biological variability.

Paper
The method implemented in this repository is detailed in the following paper:
"Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method."

Usage
To use BAMBOO batch correction, follow these steps:

Install Dependencies: Ensure that you have the required dependencies installed. These typically include Python and necessary packages like tidyverse, robustbase, openxlsx, etc.

Clone the Repository: Clone this repository to your local machine using git clone.

Run BAMBOO: Use the provided source file BAMBOO_functions to source the function to your R enviroment. Detailed usage instructions and function descriptions can be found in the source file.

Example Script: Check out the example script example.py provided in this repository to see how to apply BAMBOO batch correction to your Olink explore data.

Citation
If you find BAMBOO batch correction useful in your research, please consider citing our paper:
"Correction of batch effects in high throughput proximity extension assays for proteomic studies using bridging controls: the BAMBOO method."