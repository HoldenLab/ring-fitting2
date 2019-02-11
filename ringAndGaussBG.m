function F= ringAndGaussBG(par,imSz)

%parameters: x0, y0, r,widthRing, AmplitudeRing, bg_flat, bg_cytoplasm

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

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

R = sqrt((X-X0).^2+(Y-Y0).^2);
F = A.*exp(-((R-R0).^2./(2.*stdRing.^2)))+ bg_flat;
%defocussed gaussian cytoplasm contribution
cytoplasmBg = par(7);
%bg_cytoplasm=1300/65/2.35;%TEMP TEST
bg_cytoplasm=par(8);


F2 = cytoplasmBg.*exp(-((X-X0).^2+(Y-Y0).^2)./(2.*bg_cytoplasm.^2));

F=F+F2;
