function batchAnalyseRing_bgFitter(fileFilter,pixSz,varargin)
% function batchAnalyseRing(fileFilter,pixSz,kymoWidth,kymoBgFramespan,bleachPlotOn)
f = dir(fileFilter);
nF= numel(f);

hF= figure;
for ii = 1:nF
    fname = f(ii).name;
    fname
    try 
        ringCorrectAndFit_bgFitter(fname,pixSz,varargin{:});
    catch ME
        display(getReport(ME));
    end
end
close(hF);
%--------------------------------------------------------
function [kymo, kymoCorr] = ringCorrectAndFit_bgFitter(fname,pixSz,varargin)
%,lineProfileWidth,bleachPlotOn,,useTLFitter,frStep,doBGFG_separation)
global DEBUG_RING

ringFitArg={};
frStep = 1;
lineProfileWidth = 65;
psfFWHM = 300;
bleachPlotOn=false;
nKymoWrap=2;
nargin = numel(varargin);

ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'PsfWidthRangeNm') || strcmp(varargin{ii},'CytoplasmBG-FWHM') ...
        || strcmp(varargin{ii},'CytoplasmBG-FWHM-min') || strcmp(varargin{ii},'CytoplasmBG-FWHM-max') ...
        || strcmp(varargin{ii},'RingRadius-max') || strcmp(varargin{ii},'ZeroPadKymograph')
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+1;
    elseif strcmp(varargin{ii},'PsfFWHM')
        psfFWHM=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'LineProfleWidth')
        lineProfileWidth=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'NumKymoRepeats')
        nKymoWrap=varargin{ii+1};
        ii=ii+2;
    else
        ii=ii+1;
    end
end

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

% calculate the kymograph
[ ringStack_noBg,kymo,circFit,kymoInfo] = doBgSubAndKymo(ringStack,pixSz,lineProfileWidth,psfFWHM,frStep,ringFitArg{:});
save([savepath,filesep,fname(1:end-4),'_fitData.mat'],'circFit','kymoInfo');
%write a text file with the radius for quick reference
diamNm = round(kymoInfo(:,3))*2;
radius_fname = [savepath,filesep,fname(1:end-4),'_diamInfo','.txt'];
dlmwrite(radius_fname,diamNm);


%DEBUG CODE
%hFig =gcf;
%hold off;
%imagesc(ringStack(:,:,1))
%hold all;
%if ~useTLFitter
%    scatter(circFit.coord(:,1),circFit.coord(:,2),circFit.coord(:,4)+1);
%else
%    for ii = 1:numel(circFit);
%        scatter(circFit{ii}.coord(:,1),circFit{ii}.coord(:,2),circFit{ii}.coord(:,4)+1);
%    end
%end
%axis equal;
%colormap('gray');
%saveas(hFig, [savepath,filesep,fname(1:end-4),'_circFit.png']);
%END DEBUG

%make 2pi versions of everything for convenience
kymo_wrap = repmat(kymo,[1,nKymoWrap]);

%save everything as floats if that's what's returned
tiffwrite([savepath,filesep,fname(1:end-4),'_bgsub.tif'],ringStack_noBg);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymo.tif'],kymo);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoWrap.tif'],kymo_wrap);
