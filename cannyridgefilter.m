function [ridge ridge_nms ] = cannyridgefilter(im,blurSigma,thresh)
%function [ridge ridge_nms ] = cannyridgefilter(im,blurSigma,thresh)
%ridge filter = major eigenvalue of the image hessian

if ~exist('blurSigma','var')
    blurSigma = 1;
end
if ~exist('thresh','var')
    thresh= [1,0.5];
end



im=double(im);
%1. gaussian blur
im= imgaussfilt(im,blurSigma);

%2. ridge filter
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

%3.non-maximum suppression
[nR nC] = size(im);
ridge_nms = zeros(nR, nC);
for ii = 2:nR-1 %could make this better by padding
    for jj=2:nC-1
        if (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii+1,jj), ridge(ii-1,jj)])) || ...
            (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii,jj+1), ridge(ii,jj-1)]))
            %could do the 45 degree calculations too
            ridge_nms(ii,jj) = ridge(ii,jj);
        end
    end
end


%4. double threshold
T_high = gra
%5. Edge tracking by hysteresis
