
fname='190214_GFP-cytoplasmic.tif';

pixSz = 65; % [nm]
lineWidthNm = pixSz;
psfFWHM = 300; % [nm]

% %just test the vanilla algorithm
% batchAnalyseRing(fname,pixSz)
% display('standard alg worked ok');
%test the manual radius option
fname2 = [fname(1:end-4),'_rad300.tif'];
copyfile(fname,fname2);
batchAnalyseRing(fname2,pixSz,'SetManualRadius',300);
delete(fname2);

fname2 = [fname(1:end-4),'_rad400.tif'];
copyfile(fname,fname2);
batchAnalyseRing(fname2,pixSz,'SetManualRadius',400);
delete(fname2);
