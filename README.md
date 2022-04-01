# ExploratoryAnalysesCRinMFC

The data set and R scripts are part of a manuscript titled “Detecting Careless Responding in Multidimensional Forced-Choice Questionnaires”. 
The aim of this study was to investigate the response bias careless responding 
in the multidimensional forced-choice format. We adapted a number of indices to 
detect careless responding used in the rating scale format to the 
multidimensional forced-choice format and additionally developed several new 
indices that are unique to the multidimensional forced-choice format (The 
functions to calculate the indices are also available on GitHub [CRinMFC](https://github.com/rkupffer/CRinMFConline)). The goal of study was to investigate the performance of and correlations among these indices. Furthermore,
we investigated the extent of careless responding in MFC questionnaires as well 
as the proportion of the sample that responded carelessly.

## Coauthors
[Susanne Frick](https://github.com/susanne-frick) and
[Eunike Wetzel](https://github.com/eunike-wetzel)

## Dependencies
- [CRinMFC](https://github.com/rkupffer/CRinMFC)   
- [TirtAutomation](https://github.com/susanne-frick/TirtAutomation)   
- [MFCblockInfo](https://github.com/susanne-frick/MFCblockInfo)   

## Usage
The scripts are numbered in the order of the analyses.

1. Script 1 contains some additional information on where to access data and packages
```
1_info_access_data_packages.R
```
2. In script 2 data and syntax are prepared for Thurstonian item response theory models in Mplus
```
2_binary_coding_mplus_prep.R
```
3. Mplus analyses (either using Mplus or the package MplusAutomation)
```
3_Mplus
```
4. In script 4 the indices to detect careless responding in the multidimensional forced-choice format are calculated and saved as "4_indices.rds"
```
4_calc_indices.R
```
5. The main analyses (intercorrelations among the indices, proportion of careless responding, latent profile analyses, and analyses on the occurrence of careless responding) are part of script 5.
```
5_main_analyses.R
```
6. Additional exploratory analyses on the impact of removing careless respondents from the sample are part of script 6.
```
5_exploratory_analyses.R
```

## Licenses
Data:
[CC-BY-4.0](https://choosealicense.com/licenses/cc-by-4.0/)

Code:
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
