function [fitPar, fitIm] = fitRing(im, pixSz,psfFWHM, fixPsfFWHM,  plotOn, fixPosRad, pos, rad)
global DEBUG_RING
%parameters: x0, y0, r,width, Amplitude, BG, innerBG
RADMAX = 600/pixSz;%Like mad max

if ~exist('fixPosRad','var')
    fixPosRad=false;
end

imBlur = imgaussfilt(im,1);
%imRidge = ridgefilter(imBlur);%optional ridge filter
%im=ridgefilter(im);

%guess parameters
%segment image, ignore background
otsuThresh= otsu(imBlur);
fg = im;
fg(otsuThresh==1)=0;
%clear up any outlier pixels by only keeping the largest above thresh area
BW0=logical(fg);
BW = bwareafilt(BW0,1,'largest');
fg(BW==0)=0;


[x0C, y0C] = imCentreofMass(fg);%this sucks for uneven rings
%better use non intensity weighted avg
[x0, y0] = imBinaryCentreofMass(fg);




%r0 = 400/pixSz;
%r0 = radGuess(fg);
r0 = min(radGuess(fg),RADMAX)

d= r0*pixSz*2

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

width0 = psfFWHM/2.35/pixSz;
amplitude0 = max(fg(:));
bg0=mean(im(otsuThresh==1));
innerBG0 = 0;
cytoBg = 1300/2.35/pixSz;%manually estimated cytoplasmic gaussian contribution
cytoBgMax=2000/2.35/pixSz;
cytoBgMin=1000/2.35/pixSz;
initGuess = [x0,y0,r0,width0,amplitude0,bg0,innerBG0,cytoBg];

if ~fixPsfFWHM
    lb =[-inf, -inf,0,0,0,0,0,cytoBgMin];
    ub = [inf,inf,RADMAX,width0*4,inf,inf,inf,cytoBgMax];
else
    lb =[-inf, -inf,0,width0,0,0,0,cytoBgMin];
    ub = [inf,inf,RADMAX,width0,inf,inf,inf,cytoBgMax];
end

if fixPosRad
    lb(1:2) = pos;
    lb(3) = rad;
    ub(1:2) = pos;
    ub(3) = rad;
end
    
imSz = size(im);
%imageSizeX = imSz(2);
%imageSizeY = imSz(1);
%[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);


%DEBUG
% figure
% imagesc(im);
% pause


opt = optimoptions(@lsqcurvefit,'Display','final');

[fitPar, ~, resid] = lsqcurvefit(@(x, xdata)  ringAndGaussBG(x, xdata), ...
                         initGuess ,imSz,im, lb,ub,opt); % added resnorm 190114 kw
fitIm = ringAndGaussBG(fitPar,imSz);


ressum = sum(resid.^2);

if plotOn || (~isempty(DEBUG_RING) && DEBUG_RING==true)

    subplot(1,2,1);
    x = fitPar(1);
    y = fitPar(2);
    r = fitPar(3);
    r0*pixSz;

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
    subplot(1,2,2);
    imagesc(fitIm);
    axis equal;
    if ~isempty(DEBUG_RING) && DEBUG_RING==true
        keyboard;
    end
    pause
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
