%TEST for manual radius setting
fname = '190214_1mW_vert_1_MMStack_Pos0.ome_denoise_reg_ring2.tif'

pixSz = 65; % [nm]
psfFWHM = 300; % [nm]Check with beads

%option to fix the radius
global DEBUG_RING
DEBUG_RING=true;
manualRadius=350;
batchAnalyseRing(fname,pixSz,'SaveRawKymograph',true,'Radius',manualRadius,'CytoplasmOnlyFit',true)

figure;plot(sum(fitIm,2));
hold all;
plot(sum(im,2))
plot(sum(fitIm,1))
plot(sum(im,1))
legend('xaxis-fit','xaxis-im','yaxis-fit','yaxis-im')
