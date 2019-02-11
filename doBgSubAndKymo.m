function [ ringStack_noBg, ringKymograph, circleData, kymoInfo] = doBgSubAndKymo(ringStack,pixSzNm,lineWidthNm, psfFWHM,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings

ringStack = double(ringStack);

nFr = size(ringStack,3);

ringStack_noBg = 0.*ringStack;
for ii =1:nFr
    
    %fit each ring.
    display(['Frame: ',num2str(ii)]);
    [ringIm_noBg, ringKymoCell{ii}, circleData{ii}] = bgSubAndProfile(ringStack(:,:,ii),pixSzNm,lineWidthNm, psfFWHM,varargin{:});
    kymoSz(ii,:) = size(ringKymoCell{ii});
    rNm=circleData{ii}.r*pixSzNm;
    kymoInfo(ii,:) = [ii,rNm,kymoSz(ii,1),kymoSz(ii,2)];
    ringStack_noBg(:,:,ii) = ringIm_noBg;
end

maxKymoWidth = max(kymoSz(:,2));
ringKymograph = zeros(nFr,maxKymoWidth);
for ii = 1:nFr
    ringKymograph(ii,1:kymoSz(ii,2))=ringKymoCell{ii};
end

%------------------------------------------------------------------
function [ringIm_noBg,ringIntensity, circleData] = bgSubAndProfile(ringIm,pixSzNm,lineWidthNm, psfFWHM,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings


%fit blurred ring to the image
[fitPar,~,ringIm_noBg] = fitRing_sectored(ringIm, pixSzNm,psfFWHM, varargin{:});
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
%imagesc(ringIm);
%hold all;
%plot(Pout(:,1),Pout(:,2));
%plot(Pcirc(:,1),Pcirc(:,2));

%use this to plot a profile for all frames (lineWidth wide)
lineWidthPix = lineWidthNm/pixSzNm;
ringProfile=[];
[ringIntensity] = getProfile(ringIm_noBg,Pcirc,lineWidthPix,circ_z,circ_r,pixSzNm);
ringIntensity = ringIntensity(:)';

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


