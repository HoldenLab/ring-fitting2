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

