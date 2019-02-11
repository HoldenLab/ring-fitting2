function [fitPar fitIm bgSubIm, bgIm] = fitBlurredRing2(im, pixSz,psfFHWM, plotOn)
%parameters: x0, y0, r,width, Amplitude, BG, innerBG
if size(psfFHWM)==1
    psfFHWM_lim=[psfFHWM,psfFHWM];
else
    psfFHWM_lim=psfFHWM;
end

%guess parameters
%segment image, ignore background
T= otsu(im);
fg = im;
fg(T==1)=0;


%[x0, y0] = imCentreofMass(fg);%this sucks for uneven rings
%better use non intensity weighted avg
[x0, y0] = imBinaryCentreofMass(fg);

r0 = 400/pixSz;
widthLim = psfFHWM_lim/2.35/pixSz;
amplitude0 = max(fg(:))/2;
bg0=mean(im(T==1));
amp0_tophat =  max(fg(:));
r0_tophat = r0;
blurSz_tophat = mean(widthLim)*2;%the rings turn out really blurry (due to scattering?)
initGuess = [x0,y0,r0,mean(widthLim),amplitude0,bg0,amp0_tophat, r0_tophat,blurSz_tophat];

lb =[-inf, -inf,0,widthLim(1),0,0,0,0,0];
ub = [inf,inf,r0*4,widthLim(2),inf,inf,inf,r0*4,inf];
    
imSz = size(im);
%try initial guess with regular fit, then robust fit
fitPar_0 = lsqcurvefit(@(x, xdata)  blurredRing(x, xdata), ...
                         initGuess ,imSz,im , lb,ub);
fitPar = lsqnonlin(@(x)  blurredRing_sys(x, im),fitPar_0,lb,ub);

fitIm = blurredRing(fitPar,imSz);
%return an image with only the blurred top hat subtracted
bgPar = fitPar;
bgPar(5)=0;
bgIm = blurredRing(bgPar,imSz);
bgSubIm=im-bgIm;
ringPar = fitPar;
ringPar(7)=0;
ringIm = blurredRing(ringPar,imSz);

if plotOn
    fitPar;
    x = fitPar(1);
    y = fitPar(2);
    r = fitPar(3);

    figure; 
    imagesc(im);
    hold all;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit);
    axis equal;
    figure;
    imagesc(fitIm);
    figure;
    imagesc(bgIm);
    figure;
    imagesc(bgSubIm);
    figure;
    imagesc(ringIm);
end



%-------------------
function [x0, y0] = imCentreofMass(m)
[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

x0 = sum(vx.*x)/sum(vx);
y0 = sum(vy.*y)/sum(vy);


%------------------------------------------
function [x0, y0] = imBinaryCentreofMass(fg);
[sizey,sizex]=size(fg);
[X Y] = meshgrid(1:sizex,1:sizey);
Xring = X(fg>0);
Yring = Y(fg>0);

x0 = mean(Xring);
y0 = mean(Yring);

%-------------------------------------------
%function resid= blurredRing_sys(par,im)
function resid_w= blurredRing_sys(par,im)
%non linear function to be minimized
%parameters: x0, y0, r,width, Amplitude, BG, innerBG

 % Create a logical image of a ring with specified
% inner diameter, outer diameter center, and image size.
% First ceate the image.
imSz = size(im);
imageSizeX = imSz(2);
imageSizeY = imSz(1);
% Next create the circle in the image.
X0 = par(1);
Y0 = par(2);
R0 = par(3);
stdRing = par(4);
A = par(5);
bg = par(6);
amp_tophat=par(7);
r_tophat=par(8);
blurSz_tophat=par(9);
%the diffuse background cand be smaller diam than the ring 
%unphysical
if r_tophat<R0
    r_tophat=R0;
end
%same for blur size
if blurSz_tophat<stdRing
    blurSz_tophat=stdRing;
end

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

R = sqrt((X-X0).^2+(Y-Y0).^2);
F = A.*exp(-((R-R0).^2./(2.*stdRing.^2)))+ bg;

%now add a blurred tophat to this
tophatPar = [X0,Y0,r_tophat,blurSz_tophat,amp_tophat,0];
Ftophat = blurredTophat(tophatPar,imSz);
F = F+Ftophat;
resid = im - F;
s = mad(resid(:));% median absolute deviation of the residuals

TUNE=4.685;
r = resid./(TUNE.*s);
w = (abs(r)<1).*(1 - r.^2).^2;%BISQUARE WEIGHTING
resid_w = w.*resid; 
%imagesc(resid_w);
%pause
%-------------------------------------------
function F= blurredRing(par,imSz)

%parameters: x0, y0, r,width, Amplitude, BG, innerBG

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
bg = par(6);
amp_tophat=par(7);
r_tophat=par(8);
blurSz_tophat=par(9);
%the diffuse background cand be smaller diam than the ring 
%unphysical
if r_tophat<R0
    r_tophat=R0;
end
%same for blur size
if blurSz_tophat<stdRing
    blurSz_tophat=stdRing;
end

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

R = sqrt((X-X0).^2+(Y-Y0).^2);
F = A.*exp(-((R-R0).^2./(2.*stdRing.^2)))+ bg;

%now add a blurred tophat to this
tophatPar = [X0,Y0,r_tophat,blurSz_tophat,amp_tophat,0];
Ftophat = blurredTophat(tophatPar,imSz);
F = F+Ftophat;

%------------------------------------------------
function F= blurredTophat(par,imSz)

%parameters: x0, y0, r,edgeStdev, Amplitude, BG

 % Create a logical image of a ring with specified
% inner diameter, outer diameter center, and image size.
% First create the image.
imageSizeX = imSz(2);
imageSizeY = imSz(1);
% Next create the circle in the image.
X0 = par(1);
Y0 = par(2);
R0 = par(3);
edgeStdev = par(4);
A = par(5);
bg = par(6);

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

circlePixels = (Y - Y0).^2 ...
    + (X - X0).^2 <= R0.^2;

F = zeros(imageSizeY,imageSizeX);
F = F+ bg;
F(circlePixels) = F(circlePixels) + A;
%blur the image
F = imgaussfilt(F,edgeStdev);
%figure;imagesc(F)

