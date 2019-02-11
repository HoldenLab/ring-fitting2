function resid_w= ringAndGaussBG_robust(par,im)

%parameters: x0, y0, r,widthRing, AmplitudeRing, bg_flat,cytoplasmBg, cytoplasmBgWidth

 % Create a logical image of a ring with specified
% inner diameter, outer diameter center, and image size.
% First create the image.
imSz = size(im);
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

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

R = sqrt((X-X0).^2+(Y-Y0).^2);
F = A.*exp(-((R-R0).^2./(2.*stdRing.^2)))+ bg_flat;
%defocussed gaussian cytoplasm contribution
%cytoplasmBgWidth=1300/65/2.35;%TEMP TEST


F2 = cytoplasmBg.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*cytoplasmBgWidth.^2));

F=F+F2;

resid = im - F;
s = mad(resid(:));% median absolute deviation of the residuals

TUNE=4.685;
r = resid./(TUNE.*s);
w = (abs(r)<1).*(1 - r.^2).^2;%BISQUARE WEIGHTING
resid_w = w.*resid; 
