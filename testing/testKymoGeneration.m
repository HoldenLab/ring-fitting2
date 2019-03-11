
fname='180327_Sam1_1mW_RingHiLO_2_MMStack_Pos0.ome_denoise_reg_ring11.tif';


ringStack = double(imreadstack(fname));
pixSz = 65; % [nm]
%lineWidthNm = pixSz;
%psfFWHM = 300; % [nm]


batchAnalyseRing(fname,pixSz)
