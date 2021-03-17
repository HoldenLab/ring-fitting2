%Test analysis of VerCINI algorithm using a manually chosen radius
%This is mainly useful for analysing cells expressing cytosolic GFP 
%as a test/ calibration sample to ensure that the microscope+software give 
%even circular symettric intensity measurements around the centre of the cell
fname='190214_GFP-cytoplasmic.tif';

pixSz = 65; % [nm]

fname2 = [fname(1:end-4),'_rad300.tif'];
copyfile(fname,fname2);
verciniAnalysis(fname2,pixSz,'SetManualRadius',300);
delete(fname2);

fname2 = [fname(1:end-4),'_rad400.tif'];
copyfile(fname,fname2);
verciniAnalysis(fname2,pixSz,'SetManualRadius',400);
delete(fname2);
