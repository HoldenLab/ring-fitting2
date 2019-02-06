
fname='Denoised-reg-MAX_190128_2_MMStack_Pos0.ome-1_ring9.tif';

ringStack = double(imreadstack(fname));
pixSz = 65; % [nm]
lineWid = 10;
psfFWHM = 300; % [nm]Should be more like 300
nFrame=size(imS,3);
fRate=2;

%option to fix the radius
fixPosRad=false;
pos=[];
rad=[];
plotOn=true;
%plotOn=false
fixPsfFWHM=true;

[ringKymograph, circleData] = getRingKymo(ringStack,pixSzNm,lineWidthNm, psfFHWM)
