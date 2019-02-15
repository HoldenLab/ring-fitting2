function batchAnalyseRing(fileFilter,pixSz,varargin)
% function batchAnalyseRing(fileFilter,pixSz,varargin)
%   Performs background subtraction and kymograph fitting for vertically immobilized cells
%   Kymographs are background subtracted ie zero should equal a genuine gap in the ring
%   TODO: Optional plotting of the raw (non-bg subtracted kymograph)
% INPUTS:
%   fileFilter: Search string for files to analyse, eg '*.tif'
%   pixSz: Camera pixel size in nanometres
% OPTIONAL INPUTS:
%   'PsfFWHM', psfFWHM: Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. DEFAULT: 300 nm
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
% NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'RingRadius-max', radMax_nm: Maximum fitted ring radius. Default should hold well for WT or even most mutant Bsubtilis but change if in a different organism. If you set it too large the fitting becomes unstable for small rings. DEFAULT: 600
%   'ZeroPadKymograph', doZeroPadKymo: Add a zero row as the last row of the kymograph so that ImageJ plotting defaults to the correct contrast. DEFAULT: true
%   'FixedRadiusFit', doFixedRadiusFit: Fix the ring radius and position to the average ring position. Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. DEFAULT: true 
%   'LineProfleWidth': LineWidth: perpendicular distance over which to integrate line profile signal to improve SNR. DEFAULT:pixSz , ie 1 pixel
%   'NumKymoRepeats', nKymoWrap: Number of times to plot the kymograph side-by-side in the kymoWrap file. DEFAULT: 2
%   'SaveRawKymograph', doSaveRawKymograph: Save a non-background subtracted kymograph as well. DEFAULT, false
%
% OUTPUTS:
%   Files, in subdirectory analysed:
%       <filename>_bgsub.tif: Background subtracted ring 
%       <filename>_kymo.tif: Kymograph
%       <filename>_kymoWrap.tif: Kymograph, repeated nKymoWrap times
%       <filename>_fitData.mat: Fit results
%       <filename>_diamInfo.txt: Ring diameter, each frame

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
lineProfileWidth = pixSz;
psfFWHM = 300;
bleachPlotOn=false;
nKymoWrap=2;
doSaveRawKymograph = false;
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'PsfWidthRangeNm') || strcmp(varargin{ii},'CytoplasmBG-FWHM') ...
        || strcmp(varargin{ii},'CytoplasmBG-FWHM-min') || strcmp(varargin{ii},'CytoplasmBG-FWHM-max') ...
        || strcmp(varargin{ii},'RingRadius-max') || strcmp(varargin{ii},'ZeroPadKymograph')...
        || strcmp(varargin{ii},'FixedRadiusFit')
        ringFitArg={ringFitArg{:},varargin{ii:ii+1}};
        ii=ii+1;
    elseif strcmp(varargin{ii},'PsfFWHM')
        psfFWHM=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'LineProfileWidth')
        lineProfileWidth=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'NumKymoRepeats')
        nKymoWrap=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'SaveRawKymograph')
        doSaveRawKymograph=varargin{ii+1};
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
[ ringStack_noBg,kymo,circFit,kymoInfo, kymoraw] = doBgSubAndKymo(ringStack,pixSz,lineProfileWidth,psfFWHM,ringFitArg{:});
save([savepath,filesep,fname(1:end-4),'_fitData.mat'],'circFit','kymoInfo');
%write a text file with the radius for quick reference
diamNm = round(kymoInfo(:,3))*2;
radius_fname = [savepath,filesep,fname(1:end-4),'_diamInfo','.txt'];
dlmwrite(radius_fname,diamNm);

%make 2pi versions of everything for convenience
kymo_wrap = repmat(kymo,[1,nKymoWrap]);
kymoraw_wrap = repmat(kymoraw,[1,nKymoWrap]);

%save everything as floats if that's what's returned
tiffwrite([savepath,filesep,fname(1:end-4),'_bgsub.tif'],ringStack_noBg);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymo.tif'],kymo);
tiffwrite([savepath,filesep,fname(1:end-4),'_kymoWrap.tif'],kymo_wrap);
if doSaveRawKymograph
    tiffwrite([savepath,filesep,fname(1:end-4),'_kymoRaw.tif'],kymoraw);
    tiffwrite([savepath,filesep,fname(1:end-4),'_kymoRawWrap.tif'],kymoraw_wrap);
end

