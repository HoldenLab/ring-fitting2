function [fitPar, fitIm,ringIm_noBg] = fitRing_sectored(im, pixSz_nm,psfFWHM, varargin)
global DEBUG_RING
%parameters: x0, y0, r,width, Amplitude, BG, cytoplasmBg
NSECTOR=12;

plotOn=false;
cytoBgFWHM_nm = 1300;
cytoBgFWHMmin_nm = 1000;
cytoBgFWHMmax_nm = 2000;
radMax_nm=600;
psfWidthExtraNm= 50;%as in +/- psfFWHM
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'PsfWidthRangeNm')
        psfWidthExtraNm= varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'PlotFit')
        plotOn=true;
        ii=ii+1;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM')
        cytoBgFWHM_nm = varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM-min')
        cytoBgFWHMmin_nm = varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM-max')
        cytoBgFWHMmax_nm = varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'RingRadius-max') %best to set this pretty close to the max plausible ring radius
        radMax_nm= varargin{ii+1};
        ii=ii+2;
    else
        ii=ii+1;
    end
end

imBlur = imgaussfilt(im,1);

%guess parameters
%segment image, ignore background
otsuThresh= otsu(imBlur);
fg = im;
fg(otsuThresh==1)=0;
%clear up any outlier pixels by only keeping the largest above thresh area
BW0=logical(fg);
BW = bwareafilt(BW0,1,'largest');
fg(BW==0)=0;


%[x0C, y0C] = imCentreofMass(fg);%this sucks for uneven rings
%better use non intensity weighted avg
[x0, y0] = imBinaryCentreofMass(fg);

%%DEBUG
%hold all;
%imagesc(fg);
%th = 0:pi/50:2*pi;
%xunit = r0 * cos(th) + x0;
%yunit = r0 * sin(th) + y0;
%plot(xunit, yunit,'g');
%plot(x0, y0,'gx');
%plot(x0C, y0C,'rx');
%[x0C, y0C] = imCentreofMass(im);
%plot(x0C, y0C,'kx');
%axis equal;
%pause

imSz = size(im);

radMax = radMax_nm/pixSz_nm;%Like mad max
cytoBgWidth = cytoBgFWHM_nm/2.35/pixSz_nm;%manually estimated cytoplasmic gaussian contribution
cytoBgWidthMin=cytoBgFWHMmin_nm/2.35/pixSz_nm;
cytoBgWidthMax=cytoBgFWHMmax_nm/2.35/pixSz_nm;
width0 = psfFWHM/2.35/pixSz_nm;

r0 = min(radGuess(fg),radMax);
amplitude0 = max(fg(:));
bg0=mean(im(otsuThresh==1));
cytoplasmBg0 = max(im(:));
cytoplasmBg2_0 = mean(im(:));
cytoBgWidth2_0 =cytoBgWidth;
sectorAmp0(1:NSECTOR) = 1;

%initGuess = [x0,y0,r0,width0,amplitude0,bg0,cytoplasmBg0,cytoBgWidth,cytoplasmBg0,cytoBgWidth,sectorAmp0];
initGuess = [x0,y0,r0,width0,amplitude0,bg0,cytoplasmBg0,cytoBgWidth,cytoplasmBg2_0,cytoBgWidth2_0,sectorAmp0];

wMin = width0-psfWidthExtraNm/2.35/pixSz_nm;
wMax = width0+psfWidthExtraNm/2.35/pixSz_nm;
lb =[-inf, -inf,0,wMin,0,0,0,cytoBgWidthMin,0,cytoBgWidthMin,0.*sectorAmp0];
%ub = [inf,inf,radMax,wMax,inf,inf,inf,cytoBgWidthMax,inf,cytoBgWidthMax,1.*sectorAmp0];
ub = [inf,inf,radMax,wMax,inf,inf,inf,cytoBgWidthMax,inf,inf,1.*sectorAmp0];%the second blurred bg can be as big as you like

%imageSizeX = imSz(2);
%imageSizeY = imSz(1);
%[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);
%DEBUG
% figure
% imagesc(im);
% pause

%opt = optimoptions(@lsqcurvefit,'Display','final');
opt = optimoptions(@lsqcurvefit,'Display','off');

[fitPar] = lsqcurvefit(@(x, xdata)  ringAndGaussBG_sectored(x, xdata,NSECTOR), ...
                         initGuess ,imSz,im, lb,ub,opt); % added resnorm 190114 kw
fitIm = ringAndGaussBG_sectored(fitPar,imSz,NSECTOR);
bgPar = fitPar;
bgPar(5) = 0; %set the ring amp to 0
bgIm = ringAndGaussBG_sectored(bgPar,imSz,NSECTOR);
ringIm_noBg = im - bgIm;
ringPar = fitPar;
ringPar(7) = 0;
ringPar(9) = 0;
ringPar(6)=0;
ringIm = ringAndGaussBG_sectored(ringPar,imSz,NSECTOR);


if plotOn || (~isempty(DEBUG_RING) && DEBUG_RING==true)
    figure;

    subplot(2,2,1);
    title('raw');
    x = fitPar(1);
    y = fitPar(2);
    r = fitPar(3);
    r0*pixSz_nm;

    %figure; 
    hold off;
    %imagesc(fg);
    imagesc(im);
    colormap gray;
    hold all;
    %fit
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit,'r');
%     initial guess
    th = 0:pi/50:2*pi;
    xunit = r0 * cos(th) + x0;
    yunit = r0 * sin(th) + y0;
    plot(xunit, yunit,'g');
    axis equal;
%     legend('Fit','Initial Guess');
    %figure;
    %imagesc(fitIm);
    %axis equal;

    subplot(2,2,2);
    title('fitted image');
    imagesc(fitIm);
    colormap gray;
    axis equal;

    subplot(2,2,3);
    title('difference');
    imagesc(im-fitIm);
    colormap gray;
    axis equal;

    subplot(2,2,4);
    title('bg sub');
    imagesc(im-bgIm);
    colormap gray;
    axis equal;

        
    if ~isempty(DEBUG_RING) && DEBUG_RING==true
        figure;
        noRing=im-ringIm;
        imagesc(noRing);
        axis equal
        colormap gray
        title('Subtract the ring fit')
        tiffwrite('analysed/noRingIm.tif',noRing);
        tiffwrite('analysed/bgSub.tif',im-bgIm);
        %%REFIT THE BG
        %initGuessBg= bgPar;
        %lbBg =[fitPar(1),fitPar(2),fitPar(3),fitPar(4),0,0,0,cytoBgWidthMin,0.*sectorAmp0];
        %ubBg = [fitPar(1),fitPar(2),fitPar(3),fitPar(4),0,inf,inf,cytoBgWidthMax,0.*sectorAmp0];
        %[fitParBg] = lsqcurvefit(@(x, xdata)  ringAndGaussBG_sectored(x, xdata,NSECTOR), ...
        %                         initGuessBg ,imSz,noRing, lbBg,ubBg,opt); % added resnorm 190114 kw
        %bgIm2 =  ringAndGaussBG_sectored(fitParBg,imSz,NSECTOR);
        %figure;
        %title('Subtract the ring fit')
        %imagesc(noRing-bgIm2);
        %axis equal
        %colormap gray
        keyboard;
    else
        pause;
    end
end



%-------------------
function [x0, y0] = imCentreofMass(m)
[sizey, sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

x0 = sum(vx.*x)/sum(vx);
y0 = sum(vy.*y)/sum(vy);


% %----------------
% function Fout = blurredRingCrop(par,imSz,otsuThresh)
% 
% %parameters: x0, y0, r,width, Amplitude, BG, cytoplasmBg
% F = blurredRing(par,imSz);
% Fout = F(otsuThresh==2);

%----------------
function r0 = radGuess(fg)

[sizey,sizex]=size(fg);
[X, Y] = meshgrid(1:sizex,1:sizey);

Xring = X(fg>0);
Yring = Y(fg>0);

xRad = (max(Xring(:)) - min(Xring(:)))/2;
yRad = (max(Yring(:)) - min(Yring(:)))/2;

%check if the radii are inconsistent
if (abs(xRad-yRad)/min([xRad,yRad]))>2
    error('X and Y radius estimates are inconsistent for ring image. Fitting not reliable');
else
    r0 = mean([xRad,yRad])*0.95;
end

%------------------------------------------
function [x0, y0] = imBinaryCentreofMass(fg)
[sizey,sizex]=size(fg);
[X, Y] = meshgrid(1:sizex,1:sizey);
Xring = X(fg>0);
Yring = Y(fg>0);

x0 = mean(Xring);
y0 = mean(Yring);
