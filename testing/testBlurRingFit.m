% Author: Kevin Whitley
% Date created: 181112

% This script analyzes ring diameters from a directory.

global rawDataTheia

fname='Denoised-reg-MAX_190128_2_MMStack_Pos0.ome-1_ring9.tif';

imS = double(imreadstack(fname));
pixSz = 65; % [nm]
lineWid = 10;
psfFWHM = 300; % [nm]Check with beads
nFrame=size(imS,3);
fRate=2;

%option to fix the radius
pos=[];
rad=[];
plotOn=true;
%plotOn=false
fixPsfFWHM=true;

figure;
for jj=1:nFrame
    im=imS(:,:,jj);
    jj
    
    
    %[fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM,'PlotFit');
    %[fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM,'FitFreePsfWidth');
    %[fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM,'RingRadius-max', 2000);
    [fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM);
    rads(jj)= fitPar{jj}(3);
    stdCyto(jj) = fitPar{jj}(8);
    %pause;
end

figure
hold all
plot(rads*2*pixSz)
plot(stdCyto*2.35*pixSz);
xlabel('Frame')
ylabel('Ring diameter (nm)')
legend('ring diam','FWHM BG');
    
