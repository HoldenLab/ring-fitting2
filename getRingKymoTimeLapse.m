function [ringKymograph, circleData, kymoInfo] = getRingKymoTimeLapse(ringStack,pixSzNm,lineWidthNm, psfFHWM,nFrameStep,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings

ringStack = double(ringStack);

nFr = size(ringStack,3);

frRange = 1:nFrameStep:nFr;
if  frRange(end)<nFr
    frRange = [frRange,nFr];%make sure the entire range is covered
end

nKymoBlock = numel(frRange)-1;
for ii = 1:nKymoBlock
    ii
    frStart = frRange(ii);
    frEnd   = frRange(ii+1);
    ringSubset = ringStack(:,:,frStart:frEnd);
    [ringKymoCell{ii}, circleData{ii}] = getRingKymo(ringSubset,pixSzNm,lineWidthNm, psfFHWM,varargin{:});
    kymoSz(ii,:) = size(ringKymoCell{ii});
    rNm=circleData{ii}.r*pixSzNm;
    kymoInfo(ii,:) = [frStart,frEnd,rNm,kymoSz(ii,1),kymoSz(ii,2)];
end

maxKymoWidth = max(kymoSz(:,2));
ringKymograph = zeros(nFr,maxKymoWidth);
for ii = 1:nKymoBlock
    frStart = frRange(ii);
    frEnd   = frRange(ii+1);
    ringKymograph(frStart:frEnd,1:kymoSz(ii,2))=ringKymoCell{ii};
end
