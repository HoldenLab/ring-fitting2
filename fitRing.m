function [fitPar, fitIm,im_bgsub] = fitRing(im, pixSz_nm,psfFWHM, varargin)
% function [fitPar, fitIm,im_bgsub] = fitRing(im, pixSz_nm,psfFWHM, varargin)
%
% Fit a ring to an image.
% Blurred annulus (signal) , plus flat background plus two large width gaussians to account for diffuse and out of focus cytoplasmic background.
% Hardcoded NSECTOR (currently 12) to allow varying amplitude
% Model is:
%    X0 = par(1);
%    Y0 = par(2);
%    R0 = par(3);
%    stdRing = par(4);
%    A = par(5);
%    bg_flat = par(6);
%    cytoplasmBg = par(7);
%    cytoplasmBgWidth=par(8);
%    cytoplasmBg2 = par(9);
%    cytoplasmBgWidth2=par(10);
%    sectorAmp(1:nSector) = par(11:11+nSector-1);
%    for ii = 1:nSector
%        F_ring(sectoredImage==ii) = sectorAmp(ii)*A.*exp(-(R(sectoredImage==ii)-R0).^2./(2.*stdRing.^2));
%    end
%    %flat background contribution
%    F_bg = bg_flat;
%    %defocussed gaussian cytoplasm contribution
%    F_cyto = cytoplasmBg.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*cytoplasmBgWidth.^2));
%    %Another defocussed gaussian cytoplasm contribution
%    F_cyto2 = cytoplasmBg2.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*cytoplasmBgWidth2.^2));
%    F=F_bg+F_ring+F_cyto+F_cyto2;
%
% INPUTS:
%   im: ring Image
%   pixSz: Camera pixel size in nanometres
%   psfFWHM: Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. 
%
% OUTPUTS:
%   fitPar: As defined above.
%   fitIm: Fitted image
%   im_bgsub: Input image, im, minus the fitted background.
%       Ie im_bgsub = im - (F_bg + F_cyto + F_cyto2)
%
% OPTIONAL INPUTS:
%   'PsfFWHM', psfFWHM: Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. DEFAULT: 300 nm
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
% NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'RingRadius-max', radMax_nm: Maximum fitted ring radius. Default should hold well for WT or even most mutant Bsubtilis but change if in a different organism. If you set it too large the fitting becomes unstable for small rings. DEFAULT: 600
%   'FixedRadiusFit', fitParAvg: Fix the ring radius and position to the average ring position. Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. FITPARAVG is the result of a prior fit to an averaged ring, used to fix the positions. DEFAULT: true

global DEBUG_RING
%parameters: x0, y0, r,width, Amplitude, BG, cytoplasmBg
NSECTOR=12;

plotOn=false;
cytoBgFWHM_nm = 1300;
cytoBgFWHMmin_nm = 1000;
cytoBgFWHMmax_nm = Inf;
radMax_nm=600;
psfWidthExtraNm= 50;%as in +/- psfFWHM
doFixedRadiusFit=false;
doSetRadius=false;
doCytoOnlyFit=false;
radiusManual=NaN;
cytoFitModel = 'Gauss';
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
    elseif strcmp(varargin{ii},'FixedRadiusFit') %best to set this pretty close to the max plausible ring radius
        doFixedRadiusFit = true;
        fitParAvg= varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'Radius') %supply a manually chosen radius
        doSetRadius= true;
        radiusManual= varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoplasmOnlyFit') %force the ring amplitude to zero for cyto only fit
        doCytoOnlyFit= varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoFitModel') %force the ring amplitude to zero for cyto only fit
        % 'Gauss','Cauchy','Parametric'
        cytoFitModel= varargin{ii+1};
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

imSz = size(im);

radMax = radMax_nm/pixSz_nm;%Like mad max
width0 = psfFWHM/2.35/pixSz_nm;
if strcmp(cytoFitModel,'Gauss') %Gauss, FWHM = 2.35*w
    cytoBgWidth = cytoBgFWHM_nm/2.35/pixSz_nm;%manually estimated cytoplasmic gaussian contribution
    cytoBgWidthMin=cytoBgFWHMmin_nm/2.35/pixSz_nm;
    cytoBgWidthMax=cytoBgFWHMmax_nm/2.35/pixSz_nm;
else %cauchy distribution FWHM = 2*w
    cytoBgWidth = cytoBgFWHM_nm/2/pixSz_nm;%manually estimated cytoplasmic gaussian contribution
    cytoBgWidthMin=cytoBgFWHMmin_nm/2/pixSz_nm;
    cytoBgWidthMax=cytoBgFWHMmax_nm/2/pixSz_nm;
end

r0 = min(radGuess(fg),radMax);
amplitude0 = max(fg(:));
bg0=mean(im(otsuThresh==1));
cytoplasmBg0 = max(im(:));
cytoplasmBg2_0 = mean(im(:));
%cytoBgWidth2_0 =cytoBgWidth;
cytoBgWidth2_0 = imSz(1)/2;
sectorAmp0(1:NSECTOR) = 1;


if doFixedRadiusFit
    initGuess = fitParAvg;
    initGuess(11:11+NSECTOR-1)=ones(size(sectorAmp0));%reset the sector model
    %lb = [fitParAvg(1), fitParAvg(2), fitParAvg(3), fitParAvg(4),0,0,0,fitParAvg(8),0,fitParAvg(10),zeros(size(sectorAmp0))];
    %ub = [fitParAvg(1), fitParAvg(2), fitParAvg(3), fitParAvg(4),inf,inf,inf,fitParAvg(8),inf,fitParAvg(10),ones(size(sectorAmp0))];
    lb = [-inf,-inf, fitParAvg(3), fitParAvg(4),0,0,0,fitParAvg(8),0,fitParAvg(10),zeros(size(sectorAmp0))];
    ub = [inf,inf, fitParAvg(3), fitParAvg(4),inf,inf,inf,fitParAvg(8),inf,fitParAvg(10),ones(size(sectorAmp0))];
else
    initGuess = [x0,y0,r0,width0,amplitude0,bg0,cytoplasmBg0,cytoBgWidth,cytoplasmBg2_0,cytoBgWidth2_0,sectorAmp0];
    wMin = width0-psfWidthExtraNm/2.35/pixSz_nm;
    wMax = width0+psfWidthExtraNm/2.35/pixSz_nm;
    lb =[-inf, -inf,0,wMin,0,0,0,cytoBgWidthMin,0,cytoBgWidthMin,zeros(size(sectorAmp0))];
    ub = [inf,inf,radMax,wMax,inf,inf,inf,cytoBgWidthMax,inf,inf,ones(size(sectorAmp0))];%the second blurred bg can be as big as you like
    if doSetRadius
        radiusManualPix = radiusManual/pixSz_nm;
        initGuess(3) = radiusManualPix;
        lb(3) = radiusManualPix;
        ub(3) = radiusManualPix;
    end
    if doCytoOnlyFit
        initGuess(5) = 0;
        lb(5) = 0;
        ub(5) = 0;
    end
end

if DEBUG_RING
    opt = optimoptions(@lsqcurvefit,'Display','final');
else
    opt = optimoptions(@lsqcurvefit,'Display','off');
end

[fitPar] = lsqcurvefit(@(x, xdata)  ringAndGaussBG(x, xdata,NSECTOR,cytoFitModel), ...
                         initGuess ,imSz,im, lb,ub,opt); 

fitIm = ringAndGaussBG(fitPar,imSz,NSECTOR,cytoFitModel);
bgPar = fitPar;
bgPar(5) = 0; %set the ring amp to 0
bgIm = ringAndGaussBG(bgPar,imSz,NSECTOR,cytoFitModel);
im_bgsub = im - bgIm;
ringPar = fitPar;
ringPar(7) = 0;
ringPar(9) = 0;
ringPar(6)=0;
ringIm = ringAndGaussBG(ringPar,imSz,NSECTOR,cytoFitModel);


if plotOn || (~isempty(DEBUG_RING) && DEBUG_RING==true)
    figure;

    subplot(2,2,1);
    x = fitPar(1);
    y = fitPar(2);
    r = fitPar(3);
    r0*pixSz_nm;

    hold off;
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
    title('raw');

    subplot(2,2,2);
    imagesc(fitIm);
    colormap gray;
    axis equal;
    title('fitted image');

    subplot(2,2,3);
    imagesc(im-fitIm);
    colormap gray;
    axis equal;
    title('difference');

    subplot(2,2,4);
    imagesc(im-bgIm);
    colormap gray;
    axis equal;
    title('bg sub');

        
    if ~isempty(DEBUG_RING) && DEBUG_RING==true
        figure;
        noRing=im-ringIm;
        imagesc(noRing);
        axis equal
        colormap gray
        title('Subtract the ring fit')
        tiffwrite('analysed/noRingIm.tif',noRing);
        tiffwrite('analysed/bgSub.tif',im-bgIm);
        keyboard;
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
%------------------------------------------
function F= ringAndGaussBG(par,imSz,nSector,cytoFitModel)
%parameters: x0, y0, r,widthRing, AmplitudeRing, bg_flat, cytoplasmBgWidth

 % Create a logical image of a ring with specified
% inner diameter, outer diameter center, and image size.
% First create the image.
imageSizeX = imSz(2);
imageSizeY = imSz(1);
% Next create the circle in the image.
X0 = par(1);
Y0 = par(2);
R0 = par(3);
stdRing = par(4);
A = par(5);
bg_flat = par(6);
cytoplasmBg = par(7);
cytoplasmBgWidth=par(8);
cytoplasmBg2 = par(9);
cytoplasmBgWidth2=par(10);
sectorAmp(1:nSector) = par(11:11+nSector-1);

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

F_bg= 0.*X;
F_ring=0.*X;
F_cyto=0.*X;


R = sqrt((X-X0).^2+(Y-Y0).^2);

%sectored ring contribution
%calculate the sectored regions
%gives you an image going round in nSector sectors 1:nSector
thetaIm = atan2((Y-Y0),(X-X0));
%alternative approach
thetaLim = -pi:(2*pi)/(nSector):pi;
sectoredImage=0.*thetaIm;
for ii = 1:nSector
    sectoredImage(thetaIm>thetaLim(ii) & thetaIm<thetaLim(ii+1)) = ii;
end

for ii = 1:nSector
    F_ring(sectoredImage==ii) = sectorAmp(ii)*A.*exp(-(R(sectoredImage==ii)-R0).^2./(2.*stdRing.^2));
end
%flat background contribution
F_bg = bg_flat;
%defocussed cytoplasm contribution
if strcmp(cytoFitModel,'Gauss')
    F_cyto = cytoplasmBg.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*cytoplasmBgWidth.^2));
%DISABLE SECOND GAUSSIAN FOR NOW
%Another defocussed gaussian cytoplasm contribution
%F_cyto2 = cytoplasmBg2.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*cytoplasmBgWidth2.^2));
elseif strcmp(cytoFitModel,'Cauchy')
    %simple extension from 1d defined in kim et al 2013
    %Resolution recovery reconstruction for a Compton camera
    % using transform r^2 = x^2+y^2
    F_cyto = cytoplasmBg.*(cytoplasmBgWidth.^2./((X-X0).^2+(Y-Y0).^2+cytoplasmBgWidth.^2));

elseif strcmp(cytoFitModel,'Parametric')
end


%F=F_bg+F_ring+F_cyto+F_cyto2;
F=F_bg+F_ring+F_cyto;


