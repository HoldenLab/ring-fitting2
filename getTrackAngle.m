function [trackCirc,theta,rkymo,rtrack]= getTrackAngle(track,circleFitData)
%circleFitData is the cell of fit results (1 per frame) for each kymograph
%trackCirc = [track,theta,rtrack,rkymo];


x0=[]; 
y0=[];
r0=[];
nFr = numel(circleFitData);
for ii = 1:nFr
    x0(ii) = circleFitData{ii}.z(1);
    y0(ii) = circleFitData{ii}.z(2);
    r0(ii) = circleFitData{ii}.r;
end

nPt = size(track,1);
%trackmate zero indexes frames, matlab uses 1 index
fr = track(:,1)+1;
x = track(:,2);
y = track(:,3);

rkymo= zeros(size(x));
dx= zeros(size(x));
dy= zeros(size(x));
for ii = 1:nPt
    frCur = fr(ii);
    dx(ii,1) = x(ii) - x0(frCur);
    dy(ii,1) = y(ii) - y0(frCur);
    rkymo(ii,1) = r0(frCur);
end

theta = atan2(dy,dx);
rtrack = sqrt(dy.^2+dx.^2);

trackCirc = [track,theta,rtrack,rkymo];
%plot(fr,theta);
%pause;
