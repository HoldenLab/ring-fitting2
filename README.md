# ring-fitting2
INFORMATION 

Software for background subtraction and kymograph fitting for fluorescence microscopy image stacks of vertically immobilized bacteria using VerCINI (Vertical Cell Imaging by Nanostructured Immobilization). First published in Whitley et al, bioRxiv 2020. 
This software is provided as is and without warranty.

OVERVIEW

Performs circular kymograph analysis for VerCINI microscopy data of vertically immobilized bacteria.
Two analysis methods are provided: 
1) verciniAnalysis.m: Automated kymograph analysis where septa/ cell circumference is localized to sub-pixel precision by fitting an explicit model. The fitted background contribution is subtracted from each image, and a kymograph calculated along the circular line profile of the cell circumference. This should work for most datasets
2) manualVerciniAnalysis.m: Manual analysis where the kymograph is calculated along a circular line profile manually selected for each image stack, based on the maximum intensity projection of the stack. This is mostly used for sparse single molecule datasets if automated VerCINI analysis has failed

INSTALLATION

Add the 'ring-fitting2' directory to your MATLAB path, and save the updated path  so that it will remain installed next time you start MATLAB.

USAGE INSTRUCTIONS

To run the software, please see detailed help documentation in the main analysis functions verciniAnalysis.m, manualVerciniAnalysis.m, which can be accessed by typing 'help verciniAnalysis' on the MATLAB command line.

Test data and example analyses are included in ring-fitting2/testing


LICENSING INFORMATION

All files are distributed under the GPLv3 and (c) 2020 Seamus Holden, Newcastle University unless otherwise stated. See LICENSE.txt for full terms.

CITATION

If you use this software in work leading to a scientifc publication, please cite: 

FtsZ treadmilling is essential for Z-ring condensation and septal constriction initiation in Bacillus subtilis cell division
Kevin D. Whitley, Calum Jukes, Nicholas Tregidgo, Eleni Karinou, Pedro Almada, Ricardo Henriques, Cees Dekker, Séamus Holden
bioRxiv 2020.07.01.182006; doi: https://doi.org/10.1101/2020.07.01.182006

[![DOI](https://zenodo.org/badge/169422043.svg)](https://zenodo.org/badge/latestdoi/169422043)

