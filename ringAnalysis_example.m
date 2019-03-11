fileFilter='*.tif';
pixSz=65;
%by default, fits with fixed radius based on average image
batchAnalyseRing(fileFilter,pixSz)
%alternatively, fit with a free radius to allow for constriction
%batchAnalyseRing(fileFilter,pixSz,'FixedRadiusFit',false)
%optionally you can save the non-background subtracted kymograph
%batchAnalyseRing(fileFilter,pixSz,'SaveRawKymograph',true)
