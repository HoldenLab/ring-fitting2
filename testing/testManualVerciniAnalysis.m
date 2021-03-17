%Test analysis of manual VerCINI analysis algorithm
%Unit test and also basic introduction to the software

%1. Define the files you want to analyse
fname='180327_mNG-ftsz1.tif';
pixSzNm = 65; % [nm];

%2. run the analysis
% the file copy/delete stuff is just a hack so we dont overwrite the 
% automated analysis results. 
% you could use a similar approach if you want to try multiple different
% analysis parameters
fname2='180327_mNG-ftsz1_manual.tif';
copyfile(fname,fname2);
manualVerciniAnalysis(fname2,pixSzNm);
delete(fname2);
