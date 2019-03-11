function [trackCirc,theta,rkymo,rtrack]= getTrackAngle(track,circleFitData,pixSz)
%circleFitData is the cell of fit results (1 per frame) for each kymograph
%trackCirc = [track,theta,rtrack,rkymo];

PIXOFFSET=-1;%1 pixel shift between trackmate and matlab
circCentreOffset=0;%Matlab fitted so offset should be zero in any sane world

x0=[]; 
y0=[];
r0=[];
nFr = numel(circleFitData);
for ii = 1:nFr
    x0(ii) = circleFitData{ii}.z(1)+circCentreOffset;
    y0(ii) = circleFitData{ii}.z(2)+circCentreOffset;
    r0(ii) = circleFitData{ii}.r;
end

nPt = size(track,1);
%trackmate zero indexes frames, matlab uses 1 index
fr = track(:,1)+1;
xPix = track(:,2)/pixSz;
yPix = track(:,3)/pixSz;

rkymo= zeros(size(xPix));
dx= zeros(size(xPix));
dy= zeros(size(xPix));
for ii = 1:nPt
    frCur = fr(ii);
    dx(ii,1) = xPix(ii) - x0(frCur);
    dy(ii,1) = yPix(ii) - y0(frCur);
    rkymo(ii,1) = r0(frCur);
end

theta = atan2(dy,dx);
rtrack = sqrt(dy.^2+dx.^2);

trackCirc = [track,theta,rtrack,rkymo];
%plot(fr,theta);
%pause;
