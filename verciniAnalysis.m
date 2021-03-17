function verciniAnalysis(fileFilter,pixSz,varargin)
% Performs circular kymograph analysis for VerCINI microscopy data of vertically immobilized bacteria
% INPUTS:
%   fileFilter: Search string for files to analyse, eg '*.tif'
%   pixSz: Camera pixel size in nanometres
% OUTPUTS:
%   Files, in subdirectory analysed:
%       <filename>_bgsub.tif: Background subtracted ring 
%       <filename>_kymo.tif: Kymograph
%       <filename>_kymoWrap.tif: Kymograph, repeated nKymoWrap times
%       <filename>_fitData.mat: Fit results
%       <filename>_diamInfo.txt: Ring diameter, each frame
% OPTIONAL INPUTS (commonly used):
%   'FixedRadiusFit', true/false: Fix the ring radius and shape parameters to the average ring parameters. Note the ring centroid can still shift to allow for small drifts. 
%     Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. 
%     FITPARAVG is the result of a prior fit to an averaged ring, used to fix the positions. DEFAULT: true
%   'FixedPositionFit', true/false (default:false): Fix the ring position to the average ring position. This is important for fitting sparse data, esp single molecule, where the ring position is not well constrained. Otherwise the ring position will jump around, which is bad.  
%   'NumKymoRepeats', nKymoWrap (default:2): Number of times to plot the kymograph side-by-side in the kymoWrap file.
%   'SaveRawKymograph', true/false (default:false): Save a non-background subtracted kymograph as well. 
%   'SaveFitPNG', true/ false (default:true): Save a png of the average image overlaid with the (first) fitted circle
%   'ZeroPadKymograph', true/ false (default:true): Add a zero row as the last row of the kymograph so that ImageJ plotting defaults to the correct contrast. DEFAULT: true
% OPTIONAL INPUTS (advanced):
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:800
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:Inf
%       NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'CytoplasmOnlyFit' Dont fit a ring - useful for calibration/ testing analysis of cytoplasmic GFP samples - force the ring amplitude to zero for cyto only fit
%   'FitMaxIP', true/false (default:false): Use maximum intensity projection instead of average for the fixed radius/ position fitting
%   'LineProfileWidth', LineWidth (default:pixSz , ie 1 pixel): perpendicular distance over which to integrate line profile signal to improve SNR. 
%   'PsfFWHM', psfFWHM (default:300 nm): Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. 
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'RingRadius-max', radMax_nm (default: 600): Maximum fitted ring radius. Default should hold well for most B subtilis or E coli samples but might need to change if in a different organism. If you set it too large the fitting becomes unstable for small radius septa/ cells. 
%   'SetManualRadius', manualRadiusNm: Sets the septum radius to fixed value manualRadiusNm
%   'SetManualRadius', [manualRadiusNmMin,manualRadiusNmMax]: limits the septum radius to range defined by 2 element vector. Initial guess is average of the min and max Radius
% OPTIONAL INPUTS (developer/debug):
%   'HoughCircleGuess',true/false: (Experimental) Uses hough circle finding estimator for the initial guess. From initial testing seems more robust than the simple binary estimator. 
%
% NOTE: If the background subtration fails for some frames - slow fitting, bright bands in the kymographs - this is usually because 'FixedRadiusFit' is set to true, but the radius is changing  - try changing 'FixedRadiusFit' to false. If the radius is changign the FixedRadiusFit gives bad results as the average radius is not a good match for all frames
%
%
% EXAMPLES:
%   
%   By default, fits with fixed radius based on average image
%   >> verciniAnalysis('*.tif',pixSz)
%   Alternatively, fit with a free radius to allow for constriction
%   >> verciniAnalysis('*.tif',pixSz,'FixedRadiusFit',false)
%   Optionally you can save the non-background subtracted kymograph
%   >> verciniAnalysis('*.tif',pixSz,`'SaveRawKymograph',true)


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
doSaveFitFigure=true;
nargin = numel(varargin);
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'PsfWidthRangeNm') || strcmp(varargin{ii},'CytoplasmBG-FWHM') ...
        || strcmp(varargin{ii},'CytoplasmBG-FWHM-min') || strcmp(varargin{ii},'CytoplasmBG-FWHM-max') ...
        || strcmp(varargin{ii},'RingRadius-max') || strcmp(varargin{ii},'ZeroPadKymograph')...
        || strcmp(varargin{ii},'FixedRadiusFit') || strcmp(varargin{ii},'SetManualRadius') ...
        || strcmp(varargin{ii},'CytoplasmOnlyFit')|| strcmp(varargin{ii},'PlotFit')...
        || strcmp(varargin{ii},'ShowFitOutcome') || strcmp(varargin{ii},'FixedPositionFit') ...
        || strcmp(varargin{ii},'FitMaxIP') || strcmp(varargin{ii},'HoughCircleGuess')
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
    elseif strcmp(varargin{ii},'SaveFitPNG')
        doSaveFitFigure=varargin{ii+1};
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
Diameter_nm = round(kymoInfo(:,2))*2;
Frame= kymoInfo(:,1);
T = table(Frame,Diameter_nm);
if numel(unique(Diameter_nm))==1
    T=T(1,:);
end

radius_fname = [savepath,filesep,name,'_diamInfo','.txt'];
writetable(T,radius_fname);

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

if doSaveFitFigure
    h=figure;
    plotCircleFig(mean(ringStack,3),circFit{1});
    saveas(h,[savepath,filesep,name,'_fitResult.png']);
    close(h);
end

