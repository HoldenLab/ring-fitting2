
fname='180327_mNG-ftsz1.tif';

pixSz = 65; % [nm]
lineWidthNm = pixSz;
psfFWHM = 300; % [nm]

% %just test the vanilla algorithm
% batchAnalyseRing(fname,pixSz)
% display('standard alg worked ok');
%test the new fixPosFit
batchAnalyseRing(fname,pixSz);
display('free position fitting (default) completed ok');

fname2 = [fname(1:end-4),'_fixedPos.tif'];
copyfile(fname,fname2);
batchAnalyseRing(fname2,pixSz, 'FixedPositionFit',true)
display('fixed position fitting completed ok');
delete(fname2);
