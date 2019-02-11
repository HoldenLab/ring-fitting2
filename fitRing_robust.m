function [fitPar, fitIm] = fitRing_robust(im, pixSz_nm,psfFWHM, varargin)
global DEBUG_RING
%parameters: x0, y0, r,width, Amplitude, BG, innerBG

fixPsfFWHM=true;
plotOn=false;
cytoBgFWHM_nm = 1300;
cytoBgFWHMmin_nm = 1000;
cytoBgFWHMmax_nm = 2000;
radMax_nm=600;
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'FitFreePsfWidth')%doing this is generally a bad idea
        fixPsfFWHM=false;
        ii=ii+1;
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

radMax = radMax_nm/pixSz_nm;%Like mad max
cytoBgWidth = cytoBgFWHM_nm/2.35/pixSz_nm;%manually estimated cytoplasmic gaussian contribution
cytoBgWidthMin=cytoBgFWHMmin_nm/2.35/pixSz_nm;
cytoBgWidthMax=cytoBgFWHMmax_nm/2.35/pixSz_nm;
width0 = psfFWHM/2.35/pixSz_nm;

r0 = min(radGuess(fg),radMax);
amplitude0 = max(fg(:));
bg0=mean(im(otsuThresh==1));
innerBG0 = 0;
initGuess = [x0,y0,r0,width0,amplitude0,bg0,innerBG0,cytoBgWidth];

if ~fixPsfFWHM
    lb =[-inf, -inf,0,0,0,0,0,cytoBgWidthMin];
    ub = [inf,inf,radMax,width0*4,inf,inf,inf,cytoBgWidthMax];
else
    lb =[-inf, -inf,0,width0,0,0,0,cytoBgWidthMin];
    ub = [inf,inf,radMax,width0,inf,inf,inf,cytoBgWidthMax];
end

imSz = size(im);
%imageSizeX = imSz(2);
%imageSizeY = imSz(1);
%[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);


%DEBUG
% figure
% imagesc(im);
% pause


%opt = optimoptions(@lsqcurvefit,'Display','final');
opt = optimoptions(@lsqcurvefit,'Display','off');

%%[fitPar, ~, resid] = lsqcurvefit(@(x, xdata)  ringAndGaussBG(x, xdata), ...
%%                         initGuess ,imSz,im, lb,ub,opt); % added resnorm 190114 kw
%[fitPar, ~, resid] = lsqnonlin(@(x)  ringAndGaussBG_robust(x, im), ...
%                         initGuess ,lb,ub); 

%try an initial non-robust fit and then robust fit with fixed position parameters
[fitPar_0, ~, resid] = lsqcurvefit(@(x, xdata)  ringAndGaussBG(x, xdata), ...
                         initGuess ,imSz,im, lb,ub,opt); % added resnorm 190114 kw
initGuess2 = fitPar_0;
lb =[fitPar_0(1), fitPar_0(2),fitPar_0(3),width0,0,0,0,fitPar_0(8)];
ub =[fitPar_0(1), fitPar_0(2),fitPar_0(3),width0,inf,inf,inf,fitPar_0(8)];
[fitPar, ~, resid] = lsqnonlin(@(x)  ringAndGaussBG_robust(x, im), ...
                         initGuess2 ,lb,ub); 





fitIm = ringAndGaussBG(fitPar,imSz);
bgPar = fitPar;
bgPar(5) = 0; %set the ring amp to 0
bgIm = ringAndGaussBG(bgPar,imSz);


if plotOn || (~isempty(DEBUG_RING) && DEBUG_RING==true)

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
    axis equal;

    subplot(2,2,3);
    title('difference');
    imagesc(im-fitIm);
    axis equal;

    subplot(2,2,4);
    title('bg sub');
    imagesc(im-bgIm);
    axis equal;
    
    colormap gray;

    if ~isempty(DEBUG_RING) && DEBUG_RING==true
        keyboard;
    end
    %pause
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
% %parameters: x0, y0, r,width, Amplitude, BG, innerBG
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
