
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
%[imBgSub, ringKymograph, circleData, kymoInfo] = doBgSubAndKymo(imS,pixSz,lineWidthNm, psfFWHM)
%batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',true)
%batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',false)
batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',false,'SaveRawKymograph',true,'ShowFitOutcome',true)

load('C:\Users\nsh167\Documents\development\github\ring-fitting2\testing\analysed\180327_Sam1_1mW_RingHiLO_2_MMStack_Pos0.ome_denoise_reg_ring11_fitData.mat')
figure;hold all
w1 = fitPar(:,8)*2.35*pixSz;
w2 = fitPar(:,10)*2*pixSz;
hold all;
plot(w1);
plot(w2);
xlabel('frame');
legend('w1','w2');
saveas(gcf,'analysed/180327_Sam1_1mW_RingHiLO_1_MMStack_Pos0.ome_denoise_reg_ring1.fitted width plot.fig');

batchAnalyseRing(fname2,pixSz,'FixedRadiusFit',true,'SaveRawKymograph',true,'ShowFitOutcome',true)

