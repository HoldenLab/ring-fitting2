%Test analysis of standard VerCINI algorithm
%Unit test and also basic introduction to the software

%1. Define the files you want to analyse
%Can use wildcards to batch analyse multople files
fname='180327_*.tif'; 

%2. Define the pixel size of the image
pixSz = 65; % [nm]

%3. Run the analysis
%Results will be saved to testing/analysed
verciniAnalysis(fname,pixSz);
display('Standard VerCINI analysis algoritm completed ok');


