
fname='Denoised-1171029 Z Dynamics 2s exp 3pc 10min 150 em NH3 wf002-1_bleachCor.tif';

ringStack = double(imreadstack(fname));
pixSz = 65; % [nm]
lineWidthNm = pixSz;
psfFWHM = 300; % [nm]
nFrameStep=1

%[ringKymograph, circleData] = getRingKymo(ringStack,pixSz,lineWidthNm, psfFWHM)
[ringKymograph, circleData, kymoInfo] = getRingKymoTimeLapse(ringStack,pixSz,lineWidthNm, psfFWHM,nFrameStep);
figure;imagesc(ringKymograph);axis equal;colormap gray
