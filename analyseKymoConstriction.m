function [dt,fittable,h1] = analyseKymoConstriction(fname)
%[dt] = analyseKymoConstriction(imname);

im = imread(fname);
[h,w]=size(im);
for ii = 1:h
    [FWHM(ii), intensity(ii), widthResult] = fitProfile(im,w,ii);
end
t=1:h;



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
% d = a(3);
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
d = a(3);
wFit = piecewiseLinearWidth(a,tCrop);


dt = tConsEnd-t0;

h1=figure;
subplot(2,1,1);
hold all;
plot(FWHM,'--');
plot(tCrop,fCrop,'b-');
ylim([0 1.2*max(fCrop)])
plot(tConsEnd,0,'r*')
plot(tCrop,wFit,'r--')
ylabel('Diameter (pix)')

subplot(2,1,2);
hold all;
plot(t,intensity,'k--');
plot(tCrop,iCrop,'b-');
plot(tConsEnd,iConsEnd,'r*');
ylabel('max intensity')
xlabel('frames')

end

