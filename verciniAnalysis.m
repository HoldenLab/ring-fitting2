function verciniAnalysis(fileFilter,pixSz,varargin)
% Performs circular kymograph analysis for VerCINI microscopy data of vertically immobilized bacteria
% Septa/ cell circumference is localized to sub-pixel precision by fitting an explicit model (defined in fitRing.m) of 
% a blurred annulus two large width Gaussians to account for diffuse and out of focus cytoplasmic background, plus a baseline.
% The fitted background contribution is subtracted from each image, and a kymograph calculated along the circular line profile 
% of the cell circumference.
% By default the background subtracted image and kymograph are saved. 
% Automatically analyses all files matching "fileFilter" string
% All results are saved in 'analysed' subdirectory of the current working directory
%
% NOTE: By default, TWO copies of the kymograph are plotted side by side, ie wrapping around the cell circumference twice, 
% 0->720degrees. This is so that dynamics at the 360degree->0degree transition can still be identified
% NOTE: If the background subtraction fails for some frames - slow fitting, bright bands in the kymographs - this is usually because 'FixedRadiusFit' is set to true, but the radius is changing, or because 'FixedPositionFit' is set to true but the cell position is changing (drift) - try changing 'FixedRadiusFit' or 'FixedPositionFit' to false.
% NOTE: By default, a zero row is added as the last row of the kymograph so that ImageJ plotting defaults to the "correct" contrast for a background subtracted image, ie black = 0 counts. You can turn this off with 'ZeroPadKymograph', false
%
% INPUTS:
%   fileFilter: Search string for files to analyse, eg '*.tif'
%       Input image stacks should be cropped movies of septum/ circumferentially localized protein dynamics.
%       Input image stacks should be cropped to contain one cell only
%       Cell should not overlap image edges.
%   pixSz: Camera pixel size in nanometres
% OUTPUTS:
%   Files, in subdirectory analysed:
%       <filename>_bgsub.tif: Background subtracted vercini image stack
%       <filename>_kymo.tif: Kymograph
%       <filename>_kymoWrap.tif: Kymograph, repeated nKymoWrap times
%       <filename>_fitData.mat: Fit results
%       <filename>_diamInfo.txt: Ring diameter, each frame if constricting, otherwise only defined for first frame (same for the rest)
%
%
% OPTIONAL INPUTS (commonly used):
%   'FixedRadiusFit', true/false (default:true): Fix the ring radius and shape parameters to the average ring parameters. Note the ring centroid can still shift to allow for small drifts. 
%     Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. 
%     FITPARAVG is the result of a prior fit to an averaged ring, used to fix the positions.
%   'FixedPositionFit', true/false (default:false): Fix the ring position to the average ring position. This is set to fals by default to allow for small cell drifts during acquisition. However, for fitting sparse data, especially single molecule data, the ring position is not well constrained based only on a single image. In that case, it is important to set this option to true and use the average ring position, to avoid the ring position jumping large amounts frame to frame which is obviously bad.  This problem - fit jumping, may also occur to a lesser degree on moderately sparse non-single-molecule data, like nascent FtsZ-rings. If that is likely for your dataset, or you get weird results, try setting this option to true.
%   'NumKymoRepeats', nKymoWrap (default:2): Number of times to plot the kymograph side-by-side in the kymoWrap file.
%   'SaveRawKymograph', true/false (default:false): Save a non-background subtracted kymograph as well. 
%   'SaveFitPNG', true/ false (default:true): Save a png of the average image overlaid with the (first) fitted circle
%   'ZeroPadKymograph', true/ false (default:false): Add a zero row as the last row of the kymograph so that ImageJ plotting defaults to the correct contrast. 
% OPTIONAL INPUTS (advanced):
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:800
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:Inf
%       NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'CytoplasmOnlyFit', true/false (default:false): Dont fit a ring - useful for calibration/ testing analysis of cytoplasmic GFP samples - force the ring amplitude to zero for cyto only fit
%   'FitMaxIP', true/false (default:false): Use maximum intensity projection instead of average for the fixed radius/ position fitting. Useful for single molecule data.
%   'LineProfileWidth', LineWidth (default:pixSz , ie 1 pixel): perpendicular distance over which to integrate line profile signal to improve SNR. 
%   'PsfFWHM', psfFWHM (default:300 nm): Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. 
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'RingRadius-max', radMax_nm (default: 600): Maximum fitted ring radius. Default should hold well for most B subtilis or E coli samples but might need to change if in a different organism. If you set it too large the fitting becomes unstable for small radius septa/ cells. 
%   'SetManualRadius', manualRadiusNm: Sets the septum radius to fixed value manualRadiusNm
%   'SetManualRadius', [manualRadiusNmMin,manualRadiusNmMax]: limits the septum radius to range defined by 2 element vector. Initial guess is average of the min and max Radius
% OPTIONAL INPUTS (developer/debug):
%   'HoughCircleGuess',true/false: (Experimental) Uses hough circle finding estimator for the initial guess. From initial testing seems more robust than the simple binary estimator. 
%
%
%
% EXAMPLES:
%   
%   >> verciniAnalysis('*.tif',pixSz)
%      By default, fits with fixed radius based on average image
%
%   >> verciniAnalysis('*.tif',pixSz,'FixedRadiusFit',false)
%      Alternatively, fit with a free radius to allow for constriction
%
%   >> verciniAnalysis('*.tif',pixSz,'SaveRawKymograph',true)
%      Optionally you can save the non-background subtracted kymograph
%
%   >> verciniAnalysis('*.tif',pixSz,'FixedPositionFit', true)
%      Optionally you can fit the septum centroid position to that of the average image. Useful for single molecule or otherwise sparse data.


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

