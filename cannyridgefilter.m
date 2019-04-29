function [ridgeBW ridge nms ] = cannyridgefilter(im,varargin)
%function [ridge nms ] = cannyridgefilter(im,blurSigma,thresh)
%ridge filter = major eigenvalue of the image hessian
%based on Canny_v0.m from matlab file exchange

blurSigma = 1;
thresh= [1,0.5];
nms_con = 'horz';
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii}, 'BlurSigma')
        blurSigma=varargin{ii+1};
        ii = ii+2;
    elseif strcmp(varargin{ii}, 'Threshold')
        thresh=varargin{ii+1};
        ii = ii+2;
    elseif strcmp(varargin{ii}, 'NMS-con')
        nms_con=varargin{ii+1};
        ii = ii+2;
    else
        ii = ii+1;
    end
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

    %valley(ii) = max(e);
    ridge(ii) = -1*min(e);
end
ridge(ridge<0)=0;

%3.non-maximum suppression
[nR nC] = size(im);
nms = zeros(nR, nC);
for ii = 2:nR-1 %could make this better by padding
    for jj=2:nC-1
        if strcmp(nms_con,'horz')%this makes the most sense for kymographs
            if (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii,jj+1), ridge(ii,jj-1)])) 
                nms(ii,jj) = ridge(ii,jj);
            end
        elseif strcmp(nms_con,'4-con')
            if (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii+1,jj), ridge(ii-1,jj)])) || ...
                (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii,jj+1), ridge(ii,jj-1)])) 
                    nms(ii,jj) = ridge(ii,jj);
            end
        elseif strcmp(nms_con,'8-con')
            if (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii+1,jj+1), ridge(ii-1,jj-1)])) || ...
                (ridge(ii,jj) == max([ridge(ii,jj), ridge(ii-1,jj+1), ridge(ii+1,jj-1)]))
                    nms(ii,jj) = ridge(ii,jj);
            end
        else
            error('cannyridgefilter:unrecognizedNmsArg','Unregonized NMS connectivity argument');
        end
    end
end


%4. double threshold
Totsu = graythresh(ridge);
T_High =thresh(1) *Totsu;%or could implement absolute thresholds?
T_Low = thresh(2)*Totsu;

%5. Edge tracking by hysteresis

T_res = zeros (nR, nC);

for i = 1  : nR
    for j = 1 : nC
        if (nms(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (nms(i, j) > T_High)
            T_res(i, j) = 1;
        %Using 8-conected components
        elseif ( nms(i+1,j)>T_High || nms(i-1,j)>T_High || nms(i,j+1)>T_High || nms(i,j-1)>T_High || nms(i-1, j-1)>T_High || nms(i-1, j+1)>T_High || nms(i+1, j+1)>T_High || nms(i+1, j-1)>T_High)
            T_res(i,j) = 1;
        end;
    end;
end;

ridgeBW = logical(T_res.*255);
