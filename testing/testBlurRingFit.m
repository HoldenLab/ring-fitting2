% Author: Kevin Whitley
% Date created: 181112

% This script analyzes ring diameters from a directory.

global rawDataTheia

fname='Denoised-reg-MAX_190128_2_MMStack_Pos0.ome-1_ring9.tif';

imS = double(imreadstack(fname));
pixSz = 65; % [nm]
lineWid = 10;
psfFWHM = 300; % [nm]Should be more like 300
nFrame=size(imS,3);
fRate=2;

%option to fix the radius
fixPosRad=false;
pos=[];
rad=[];
plotOn=true;
%plotOn=false;
fixPsfFWHM=true;

figure;
for jj=1:nFrame
    im=imS(:,:,jj);
    jj
    
    
    [fitPar{jj}, fitIm] = fitRing(im, pixSz,psfFWHM, fixPsfFWHM,  plotOn, fixPosRad, pos, rad);
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
    
