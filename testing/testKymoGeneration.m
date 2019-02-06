
fname='Denoised-ftsz sam1  hilo4593 3pc 1s em150 001.reg-1.tif';

ringStack = double(imreadstack(fname));
pixSz = 65; % [nm]
lineWidthNm = pixSz;
psfFWHM = 300; % [nm]


[ringKymograph, circleData] = getRingKymo(ringStack,pixSz,lineWidthNm, psfFWHM)
%nFrame=size(ringStack,3);
%fRate=2;
figure;
imagesc(ringKymograph)
axis equal;
colormap gray
