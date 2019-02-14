
global DEBUG_RING
%DEBUG_RING=true;
DEBUG_RING=false;

fname1='180327_Sam1_1mW_RingHiLO_1_MMStack_Pos0.ome_denoise_reg_ring1.tif';
fname2='180327_Sam1_1mW_RingHiLO_2_MMStack_Pos0.ome_denoise_reg_ring11.tif';

pixSz = 65; % [nm]
psfFWHM = 300; % [nm]Check with beads
fRate=2;

%option to fix the radius
pos=[];
rad=[];
plotOn=true;
%plotOn=false
lineWidthNm = pixSz;
%[imBgSub, ringKymograph, circleData, kymoInfo] = doBgSubAndKymo(imS,pixSz,lineWidthNm, psfFWHM)
batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',true)
%batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',false)
batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',true,'SaveRawKymograph',true)
