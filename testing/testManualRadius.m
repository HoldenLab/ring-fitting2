%TEST for manual radius setting
fname = '190214_1mW_vert_1_MMStack_Pos0.ome_denoise_reg_ring2.tif'

pixSz = 65; % [nm]
psfFWHM = 300; % [nm]Check with beads

%option to fix the radius
manualRadius=350;
batchAnalyseRing(fname,pixSz,'SaveRawKymograph',true,'Radius',manualRadius)%I dont think bg subtracted data will necessarily make lots of sense with the cyto gfp data?

