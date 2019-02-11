
global DEBUG_RING
DEBUG_RING=true;

fname1='ring_hilo\180327_Sam1_1mW_RingHiLO_1_MMStack_Pos0.ome_denoise_reg_ring1.tif';
fname2='ring_hilo\180327_Sam1_1mW_RingHiLO_2_MMStack_Pos0.ome_denoise_reg_ring11.tif';

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
fixPsfFWHM=true;

jj=1
im=imS(:,:,jj);
%[fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM);
%rads(jj)= fitPar{jj}(3);
%stdCyto(jj) = fitPar{jj}(8);
%[fitPar{jj}, fitIm] = fitRing_robust(im, pixSz,psfFWHM);
%rads(jj)= fitPar{jj}(3);
%stdCyto(jj) = fitPar{jj}(8);
[fitPar{jj}, fitIm] = fitRing_sectored(im, pixSz,psfFWHM,'PsfWidthRangeNm',1000);
rads(jj)= fitPar{jj}(3);
stdCyto(jj) = fitPar{jj}(8);

