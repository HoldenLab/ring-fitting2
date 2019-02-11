
global DEBUG_RING
DEBUG_RING=false;

fname1='180327_Sam1_1mW_RingHiLO_1_MMStack_Pos0.ome_denoise_reg_ring1.tif';
fname2='180327_Sam1_1mW_RingHiLO_2_MMStack_Pos0.ome_denoise_reg_ring11.tif';

imS = double(imreadstack(fname1));
pixSz = 65; % [nm]
psfFWHM = 300; % [nm]Check with beads
nFrame=size(imS,3);
fRate=2;

%option to fix the radius
pos=[];
rad=[];
plotOn=true;
%plotOn=false
lineWidthNm = pixSz;
%[imBgSub, ringKymograph, circleData, kymoInfo] = doBgSubAndKymo(imS,pixSz,lineWidthNm, psfFWHM)
batchAnalyseRing_bgFitter(fname2,pixSz)
