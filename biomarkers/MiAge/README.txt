"main_R_script.r" : This is the main R code which calculates mitotic ages of input samples according to our method. Only this code needs to be run.

The input for this code, "input_HM450_methyl_data.txt" is a matrix of HM450 methylation values with each row representing each CpG site and each column representing each sample whose mitotic ages will be estimated. The row names of CpG sites and column names corresponding to sample names must be given.                                             

"function_library.r" : Library of all functions used in the calculation of pvalues
"site_specific_parameters.Rdata": This Rdata contains estimates for the site-specific parameters of the epigenetic clock sites
"Additional_File1.csv" : information for the epigenetic clock sites





