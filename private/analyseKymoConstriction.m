function [dt,fittable,h1,diam0] = analyseKymoConstriction(fname,pixsz)
%[dt] = analyseKymoConstriction(imname);

if ~exist('pixsz')
    pixsz=65;
end

MAXDIAM = 1000;
MINDIAM = 50;
im = imread(fname);

[nfr,w]=size(im);

immin=1;
immax=nfr;
cropname = [fname(1:end-4),'_frlim.txt'];
if exist(cropname,'file')
    croplim=importdata(cropname);
    immin = croplim(1);
    if numel(croplim)>1
        immax= croplim(2);
    end
else
    %make the file
    frlim=[immin;immax];
    dlmwrite(cropname,frlim)
end

kk=1;
for ii = immin:immax
    [FWHM(kk), intensity(kk), widthResult] = fitProfile(im,w,ii);
    kk=kk+1;
end
FWHM=FWHM.*pixsz;
t=immin:immax;
badIdx = FWHM>MAXDIAM | FWHM<MINDIAM;
FWHM(badIdx)=[];
intensity(badIdx)=[];
t(badIdx)=[];

% need a way to delete all the crap after constriction finishes
% Constriction end is just after intensity peak
[iMax, idxMax]= max(intensity);
tMax= t(idxMax);

[iMin, idx]= min(intensity);
tMin=t(idx);
iConsEnd = (iMax+iMin)/2;
tCrop= t(idxMax:end);
iCrop = intensity(idxMax:end);

iDiff = iCrop-iConsEnd;
iConsEnd = iCrop(find(iDiff<0,1)-1);
tConsEnd = tCrop(find(iDiff<0,1)-1);

%format output as table
fittable = table(t',FWHM',intensity','VariableNames',{'Frame','FWHM','maxintensity'});

%fit the width function to the data
% t0 = a(1);
% t1 = a(2);
% diam0 = a(3);
tCrop= t(1:idxMax);
fCrop = FWHM(1:idxMax);
iCrop = intensity(1:idxMax);
a0= [min(tCrop),tConsEnd,max(fCrop)];
lb=[0,0,0];
ub=[inf,inf,inf];
options = optimoptions('lsqcurvefit','Display','off');
a = lsqcurvefit(@piecewiseLinearWidth,a0,tCrop,fCrop,lb,ub,options);
t0 = a(1);
t1 = a(2);
diam0 = a(3);
wFit = piecewiseLinearWidth(a,tCrop);

dt = tConsEnd-t0;

h1=figure;
subplot(2,1,1);
hold all;
plot(t,FWHM,'--');
plot(tCrop,fCrop,'b-');
ylim([0 1.2*max(fCrop)])
plot(tConsEnd,0,'r*')
plot(tCrop,wFit,'r--')
ylabel('Diameter (pix)')
xlim([1,nfr]);

subplot(2,1,2);
hold all;
plot(t,intensity,'k--');
plot(tCrop,iCrop,'b-');
plot(tConsEnd,iConsEnd,'r*');
ylabel('max intensity')
xlabel('frames')
xlim([1,nfr]);

end

