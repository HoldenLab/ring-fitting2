fileFilter='*.tif';
pixelSize=64;
%by default, fits with fixed radius based on average image
batchAnalyseRing(fileFilter,pixSz)
%alternatively, fit with a free radius to allow for constriction
%batchAnalyseRing(fileFilter,pixSz,'FixedRadiusFit',true)
