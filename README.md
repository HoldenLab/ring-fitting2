# ring-fitting2
INFORMATION 
MATLAB software for background subtraction and kymograph fitting for fluorescence microscopy of vertically immobilized bacteria using VerCINI (Vertical Cell Imaging by Nanostructured Immobilization). Associated with Whitley et al, bioRxiv 2020.  
This code is provided as is and without warranty to support reproducibility and open science. It is in-house code, minimally documented, warts and all. We aim increase the user friendliness and documentation of this software in the near future - watch this space!

INSTRUCTIONS
batchAnalyseRing.m : USE THIS to run everything. Analyses one or more ring movies to background subtract and plot kymogrpaphs. See documentation in the file.   
ringAnalysis_example.m: Example script showing how to use this fitting library.
doBgSubAndKymo.m: Background subtract and plot kymogrpaphs for a single ring movie. See documentation in the file.
fitRing.m: Fit a single ring image to a blurred ring + background model. See documentation in the file.              
testing\ :  Various VerCINI data examples and test scripts.

INSTALLATION
Add the 'ring-fitting2' directory to your MATLAB path, and save the updated path  so that it will remain installed next time you start MATLAB.

LICENSING INFORMATION
All files are distributed under the GPLv3 and (c) 2020 Seamus Holden, Newcastle University unless otherwise stated. See LICENSE.txt for full terms.
