function [ringKymograph, circleData, kymoInfo] = getRingKymoTimeLapse(ringStack,pixSzNm,lineWidthNm, psfFHWM,nFrameStep,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings

ringStack = double(ringStack);

nFr = size(ringStack,3);

frStart= 1:nFrameStep:nFr;
nKymoBlock = numel(frStart);
frEnd   = frStart+nFrameStep-1;
if frEnd(end)>nFr
    frEnd(end)=nFr;
end
keyboard
for ii = 1:nKymoBlock
    display(['Frame: ',num2str(frStart(ii))]);
    ringSubset = ringStack(:,:,frStart(ii):frEnd(ii));
    [ringKymoCell{ii}, circleData{ii}] = getRingKymo(ringSubset,pixSzNm,lineWidthNm, psfFHWM,varargin{:});
    kymoSz(ii,:) = size(ringKymoCell{ii});
    rNm=circleData{ii}.r*pixSzNm;
    kymoInfo(ii,:) = [frStart(ii),frEnd(ii),rNm,kymoSz(ii,1),kymoSz(ii,2)];
end

maxKymoWidth = max(kymoSz(:,2));
ringKymograph = zeros(nFr,maxKymoWidth);
for ii = 1:nKymoBlock
    ringKymograph(frStart(ii):frEnd(ii),1:kymoSz(ii,2))=ringKymoCell{ii};
end
