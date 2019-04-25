function [ridge valley gxx gyy ] = ridgefilter(im,blurSigma)
%ridge filter = major eigenvalue of the image hessian
if ~exist('blurSigma','var')
    blurSigma = 1;
end
im= imgaussfilt(im,blurSigma);

[gx, gy] = gradient(double(im));
[gxx, gxy] = gradient(gx);
[gxy, gyy] = gradient(gy);

valley = zeros(size(im));
ridge = zeros(size(im));
for ii = 1:numel(im)
    A = [gxx(ii), gxy(ii);...
         gxy(ii), gyy(ii)];

    e = eig(A);

    valley(ii) = max(e);
    ridge(ii) = -1*min(e);
end
