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

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

R = sqrt((X-X0).^2+(Y-Y0).^2);
F = A.*exp(-((R-R0).^2./(2.*stdRing.^2)))+ bg;

%now add a blurred tophat to this
if numel(par)>6
    innerBg = par(7);
    tophatPar = [X0,Y0,R0,stdRing,innerBg,0];
    Ftophat = blurredTophat(tophatPar,imSz);
    F = F+Ftophat;
end
