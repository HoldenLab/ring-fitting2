function [BW ridge] = ridgefinder(im, blurSigma)
if ~exist('blurSigma','var')
    blurSigma = 2;
end

% blur
if blurSigma> 0
    imBlur = imgaussfilt(im,blurSigma);
else
    imBlur = im;
end

ridge = ridgefilter(imBlur);

T = adaptthresh(ridge);
BW = imbinarize(ridge,T);

