function [ringKymograph, circleData] = getRingKymoWide2(ringStack,pixSzNm,lineWidthNm, psfFHWM)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings

ringStack = double(ringStack);
ringStack_avg = mean(ringStack,3);

%fit blurred ring to the image
fixPsfFWHM = true;
plotOn = true;
fitPar= fitBlurRingRidge(ringStack_avg,pixSzNm,psfFHWM,fixPsfFWHM,plotOn);
x=fitPar(1);
y=fitPar(2);
circ_z=[x,y];
circ_r=fitPar(3);
% calculate a single pixel of circumference and use that as the 
% spacing to sample the circle
theta_step = 1/circ_r;
% non clockwise non-top thetat really confusing
%define clockwise from twelve-o-clock
theta = 0:theta_step:2*pi;
%theta = pi/2:-theta_step:-3/2*pi;
distStep =cumsum([0,ones(1,numel(theta)-1)*abs(theta(2)-theta(1))]);
Pcirc = [circ_z(1)+circ_r*cos(theta(:)), circ_z(2)+circ_r*sin(theta(:)), theta(:), distStep(:)];

%%DEBUG
%figure;
%imagesc(ringStack_avg);
%hold all;
%plot(Pout(:,1),Pout(:,2));
%plot(Pcirc(:,1),Pcirc(:,2));

%use this to plot a profile for all frames (lineWidth wide)
lineWidthPix = lineWidthNm/pixSzNm;
Options=struct;
ringProfile=[];
nFr = size(ringStack,3);
for ii = 1:nFr
    I = ringStack(:,:,ii);
    [ringIntensity] = getProfile(I,Pcirc,lineWidthPix,circ_z,circ_r,pixSzNm);
    ringKymograph(ii,:) = ringIntensity(:)';
end

circleData.coord = Pcirc;
circleData.z = circ_z;
circleData.r = circ_r;

%------------------------------------------------------------------
function [profileIntensityWide, profileIntensity1pix] = getProfile(I,samplePts,lineWidthPix,circ_z,circ_r,pixSzNm);

x = samplePts(:,1);
y= samplePts(:,2);
nPt = numel(x);

STEPSZ =0.1;
%nLineStep is in each direction +/-
nLineStep = round(lineWidthPix/2/STEPSZ);

%for each point, get the intensity
[Xim,Yim] = meshgrid(1:size(I,2),1:size(I,1));
for ii = 1:nPt
    %get a perpendicular line defining the width to integrate along
    XC = [x(ii),y(ii)];
    dX = [XC(1) - circ_z(1),XC(2) - circ_z(2)];
    dXnorm = dX./sqrt(sum(dX.^2));%normalize
    dL = -nLineStep*.1:.1:nLineStep*.1;
    Xline = [repmat(XC',1,numel(dL)) + dXnorm'*dL]';
    %%interpolate the intensity to sub-pixel precision
    lineProfile = interp2(Xim,Yim,I,Xline(:,1),Xline(:,2),'cubic');
    profileIntensityWide(ii) =trapz(dL,lineProfile);
    profileIntensity1pix(ii) = interp2(Xim,Yim,I,x(ii),y(ii),'cubic');
    %figure(1);
    %hold off
    %imagesc(I);
    %hold all;
    %plot(XC(1),XC(2),'x');
    %plot(Xline(:,1),Xline(:,2),'-');
    %figure(2);
    %plot(dL,lineProfile);
    %pause
end

