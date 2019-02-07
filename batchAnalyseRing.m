function batchAnalyseRing(fileFilter,pixSz,varargin)
% function batchAnalyseRing(fileFilter,pixSz,kymoWidth,kymoBgFramespan,bleachPlotOn)
f = dir(fileFilter);
nF= numel(f);
for ii = 1:nF
    fname = f(ii).name;
    fname
    try 
        ringCorrectAndFit(fname,pixSz,varargin{:});
    catch ME
        display(getReport(ME));
    end
end
%--------------------------------------------------------
function [kymoRaw, kymoCorr] = ringCorrectAndFit(fname,pixSz,varargin)
%,lineProfileWidth,bleachPlotOn,,useTLFitter,frStep,doBGFG_separation)
global DEBUG_RING

ringFitArg={};
useTLFitter = false;
frStep = 1;
lineProfileWidth = 65;
psfFWHM = 300;
bleachPlotOn=false;
nargin = numel(varargin);

ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'FitFreePsfWidth')%doing this is generally a bad idea
        ringFitArg={ringFitArg{:},varargin{ii}};
        ii=ii+1;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM')
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM-min')
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+2;
    elseif strcmp(varargin{ii},'CytoplasmBG-FWHM-max')
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+2;
    elseif strcmp(varargin{ii},'RingRadius-max') %best to set this pretty close to the max plausible ring radius
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+2;
    elseif strcmp(varargin{ii},'TimeLapseFitter')
        useTLFitter = true;
        frStep = varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'PsfFWHM')
        psfFWHM=varargin{ii+1};
        ii=ii+1;
    elseif strcmp(varargin{ii},'LineProfleWidth')
        lineProfileWidth=varargin{ii+1};
        ii=ii+1;
    else
        ii=ii+1;
    end
end

%if ~exist('bleachPlotOn','var')
%    bleachPlotOn=false;
%end
%if ~exist('useTLFitter','var')
%    useTLFitter=false;
%    frStep=[];
%end

[pathstr,name,ext]=fileparts(fname);
if isempty(pathstr)
    pathstr='.';%matlab looks on all pths not just current path if dont do this
end

if ~isempty(DEBUG_RING) && DEBUG_RING==true
    savepath = fullfile(pathstr,'test_analysis');
else
    savepath = fullfile(pathstr,'analysed');
end

if ~exist(savepath,'dir')
    mkdir(savepath);
end

ringStack = imreadstack(fname);

[ringStackCorr,  ringStackData_bgSubOnly, tau_fg, tau_bg,baseLineBG] = otsuBleachCorrect(ringStack,bleachPlotOn);
bleachFit.tau_fg= tau_fg;
bleachFit.tau_bg= tau_bg;
bleachFit.baseLineBG= baseLineBG;

% calculate the kymograph
if useTLFitter
    %[kymoRaw,circFitRaw,kymoInfoRaw] = getRingKymoTringStackeLapse(ringStack,pixSz,lineProfileWidth,psfFWHM,frStep);
    %[kymoCorr, circFitCorr, kymoInfoCorr] = getRingKymoTringStackeLapse(ringStackCorr,pixSz,lineProfileWidth,psfFWHM,frStep);
    %save([savepath,filesep,fname(1:end-4),'_fitData.mat'],'bleachFit','circFitCorr','circFitRaw','kymoInfoRaw','kymoInfoCorr');
    %%write a blank text file with the radius for quick reference
    %diamNm = round(kymoInfoCorr(:,3))*2;
    %radius_fname = [savepath,filesep,fname(1:end-4),'_diamInfo','.txt'];
    %dlmwrite(radius_fname,diamNm);
else
    [kymoRaw] = getRingKymo(ringStack,pixSz,lineProfileWidth,psfFWHM,ringFitArg{:});
    [kymoCorr circFit] = getRingKymo(ringStackCorr,pixSz,lineProfileWidth,psfFWHM,ringFitArg{:});
    save([savepath,filesep,fname(1:end-4),'_fitData.mat'],'bleachFit','circFit');
    %write a blank text file with the radius for quick reference
    f = fopen([savepath,filesep,fname(1:end-4),'_diam',num2str(round(2*circFit.r*pixSz)),'.txt'],'w');
    fclose(f);
end

%figure;
hFig =gcf;
hold off;
imagesc(ringStackCorr(:,:,1))
hold all;
if ~useTLFitter
    scatter(circFit.coord(:,1),circFit.coord(:,2),circFit.coord(:,4)+1);
else
    for ii = 1:numel(circFitCorr);
        circFit = circFitCorr{ii};
        scatter(circFit.coord(:,1),circFit.coord(:,2),circFit.coord(:,4)+1);
    end
end
axis equal;
colormap('gray');
saveas(hFig, [savepath,filesep,fname(1:end-4),'_circFit.png']);


%make 2pi versions of everything for convenience
kymoRaw_wrap = repmat(kymoRaw,[1,2]);
kymoCorr_wrap = repmat(kymoCorr,[1,2]);

%save everything as floats if that's what's returned
tiffwrite([savepath,filesep,fname(1:end-4),'_bleachCor.tif'],ringStackCorr);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoRaw.tif'],kymoRaw);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoRaw_wrap.tif'],kymoRaw_wrap);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoCorr.tif'],kymoCorr);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoCorr_wrap.tif'],kymoCorr_wrap);
