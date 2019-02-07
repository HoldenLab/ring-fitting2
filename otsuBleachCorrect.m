function [imDataCorr, imData_bgSubOnly, tau_fg, tau_bg,baseLineBG, curveFG, curveBG, hFig] = otsuBleachCorrect(imData,plotOn,nOtsuLevel)
%corrects background and corrects foreground bleaching 
% uses otsus method to separate signal and background
% returns double
imData = double(imData);

if ~exist('plotOn','var')
    plotOn = false;
end
if~exist('nOtsuLevel','var')
    nOtsuLevel=2;
end

%For each frame, 
    %do Otsu
    %estimate FG, BG (mean)
nFrame = size(imData,3);

for ii = 1:nFrame
    a = imData(:,:,ii);
    T= otsu(a,nOtsuLevel);
    fg(ii) = double(mean(a(T>=2)));
    bg(ii) = double(mean(a(T==1)));
end
t = 1:nFrame;

ft = fittype('a*exp(-b*x)+c');
[curveBG,gof1] = fit(t',bg',ft,'StartPoint',[max(bg), 1, 0],'Lower',[0 0 0]);

% subtract the fitted BG
bgFit = curveBG(t)';
fg_z = fg - bgFit;

%try fitting offset corrected FG using exp with no offset
[curveFG,gofFG] = fit(t',fg_z','exp1','Lower',[0 -inf], 'Upper',[inf, 0]);

tau_bg = curveBG.b^-1;
baseLineBG = curveBG.c;
tau_fg = -curveFG.b^-1;
if plotOn
    hFig = gcf;
    hold off;
    plot(t,fg,'o');
    hold all;
    plot(t,bg,'o');
    plot(t,fg_z,'o');
    plot(t,curveBG(t));
    plot(t,curveFG(t));
    bGtext = ['Fitted BG, tau=',num2str(round(tau_bg)),'fr, offset ',num2str(round(baseLineBG))];
    fGtext = ['Fitted FG, tau=',num2str(round(tau_fg)),' fr'];
    legend('Foreground','Background','FG - BG', bGtext,fGtext);
    %with offset - whats the tau look like?
    %[curveFG2,gofFG2] = fit(t',fg_z',ft,'StartPoint',[max(fg_z), 1, 0])
    %plot(t,curveFG2(t));
    %tau_fg2 = curveFG2.b^-1
    %imData(:,:,1);
    %T= otsu(a,nOtsuLevel);
    %figure;imagesc(T);
    %title('frame 1 otsu')
else
    hFig=[];
end

%rescale the signal to the unbleached value
fg_fit = curveFG(t);
rescaleVal = fg_fit(1)./fg_fit;

% Bleach corr FG
imDataCorr = imData;
for ii = 1:nFrame
    imDataCorr(:,:,ii) = imDataCorr(:,:,ii) - bgFit(ii);%remove background
    imData_bgSubOnly(:,:,ii) = imDataCorr(:,:,ii);

    imDataCorr(:,:,ii) = imDataCorr(:,:,ii)*rescaleVal(ii); %rescale signal
end


