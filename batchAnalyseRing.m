function batchAnalyseRing(fileFilter,pixSz,varargin)
% function batchAnalyseRing(fileFilter,pixSz,varargin)
%   Performs background subtraction and kymograph fitting for vertically immobilized cells
%   Kymographs are background subtracted ie zero should equal a genuine gap in the ring
% INPUTS:
%   fileFilter: Search string for files to analyse, eg '*.tif'
%   pixSz: Camera pixel size in nanometres
% OPTIONAL INPUTS:
%   'PsfFWHM', psfFWHM: Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. DEFAULT: 300 nm
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:800
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:Inf
%   'CytoplasmOnlyFit' Dont fit a ring - useful for cytoGFP - force the ring amplitude to zero for cyto only fit
% NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'RingRadius-max', radMax_nm: Maximum fitted ring radius. Default should hold well for WT or even most mutant Bsubtilis but change if in a different organism. If you set it too large the fitting becomes unstable for small rings. DEFAULT: 600
%   'ZeroPadKymograph', doZeroPadKymo: Add a zero row as the last row of the kymograph so that ImageJ plotting defaults to the correct contrast. DEFAULT: true
%   'FixedRadiusFit', true/false: Fix the ring radius and shape parameters to the average ring parameters. Note the ring centroid can still shift to allow for small drifts. 
%     Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. 
%     FITPARAVG is the result of a prior fit to an averaged ring, used to fix the positions. DEFAULT: true
%   'FixedPositionFit', true/false: Fix the ring position to the average ring position. This is important for fitting sparse data, esp single molecule, where the ring position is not well constrained.
%     otherwise the ring position will jump around, which is bad.  DEFAULT: false 
%   'LineProfileWidth': LineWidth: perpendicular distance over which to integrate line profile signal to improve SNR. DEFAULT:pixSz , ie 1 pixel
%   'NumKymoRepeats', nKymoWrap: Number of times to plot the kymograph side-by-side in the kymoWrap file. DEFAULT: 2
%   'SaveRawKymograph', true/false: Save a non-background subtracted kymograph as well. DEFAULT, false
%
% NOTE: If the background subtration fails for some frames - slow fitting, bright bands in the kymographs - this is usually because 'FixedRadiusFit' is set to true, but the radius is changing  - try changing 'FixedRadiusFit' to false. If the radius is changign the FixedRadiusFit gives bad results as the average radius is not a good match for all frames
%
% OUTPUTS:
%   Files, in subdirectory analysed:
%       <filename>_bgsub.tif: Background subtracted ring 
%       <filename>_kymo.tif: Kymograph
%       <filename>_kymoWrap.tif: Kymograph, repeated nKymoWrap times
%       <filename>_fitData.mat: Fit results
%       <filename>_diamInfo.txt: Ring diameter, each frame
%
% EXAMPLES:
%   
%   By default, fits with fixed radius based on average image
%   >> batchAnalyseRing('*.tif',pixSz)
%   Alternatively, fit with a free radius to allow for constriction
%   >> batchAnalyseRing('*.tif',pixSz,'FixedRadiusFit',false)
%   Optionally you can save the non-background subtracted kymograph
%   >> batchAnalyseRing('*.tif',pixSz,`'SaveRawKymograph',true)
%TODO fix the diameter file bug

f = dir(fileFilter);
nF= numel(f);

hF= figure;
for ii = 1:nF
    file = f(ii).name;
    folder=f(ii).folder;
    fname = fullfile(folder,file)
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
doSaveRawKymograph = true;
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'PsfWidthRangeNm') || strcmp(varargin{ii},'CytoplasmBG-FWHM') ...
        || strcmp(varargin{ii},'CytoplasmBG-FWHM-min') || strcmp(varargin{ii},'CytoplasmBG-FWHM-max') ...
        || strcmp(varargin{ii},'RingRadius-max') || strcmp(varargin{ii},'ZeroPadKymograph')...
        || strcmp(varargin{ii},'FixedRadiusFit') || strcmp(varargin{ii},'Radius') ...
        || strcmp(varargin{ii},'CytoplasmOnlyFit')|| strcmp(varargin{ii},'PlotFit')...
        || strcmp(varargin{ii},'ShowFitOutcome') || strcmp(varargin{ii},'FixedPositionFit')
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
[ ringStack_noBg,kymo,circFit,kymoInfo, kymoraw,fitPar] = doBgSubAndKymo(ringStack,pixSz,lineProfileWidth,psfFWHM,ringFitArg{:});
save([savepath,filesep,name,'_fitData.mat'],'circFit','kymoInfo','fitPar');
%write a text file with the radius for quick reference
diamNm = round(kymoInfo(:,3))*2;
radius_fname = [savepath,filesep,name,'_diamInfo','.txt'];
dlmwrite(radius_fname,diamNm);

%make 2pi versions of everything for convenience
kymo_wrap = repmat(kymo,[1,nKymoWrap]);
kymoraw_wrap = repmat(kymoraw,[1,nKymoWrap]);

%save everything as floats if that's what's returned
tiffwrite([savepath,filesep,name,'_bgsub.tif'],ringStack_noBg);
tiffwrite([savepath,filesep,name,'_kymo.tif'],kymo);
tiffwrite([savepath,filesep,name,'_kymoWrap.tif'],kymo_wrap);
if doSaveRawKymograph
    tiffwrite([savepath,filesep,name,'_kymoRaw.tif'],kymoraw);
    tiffwrite([savepath,filesep,name,'_kymoRawWrap.tif'],kymoraw_wrap);
end

