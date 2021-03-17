function manualVerciniAnalysis(fileFilter,pixSz,nKymoWrap)
% Performs manual circular kymograph analysis for VerCINI microscopy data of vertically immobilized bacteria
% Kymograph is calculated along a circular line profile manually selected for each image stack, based on the maximum
% intensity projection of the stack
% This is mostly used for sparse single molecule datasets if automated VerCINI analysis has failed
% Automatically analyses all files matching "fileFilter" string
%
% NOTE: By default, TWO copies of the kymograph are plotted side by side, ie wrapping around the cell circumference twice, 
% 0->720degrees. This is so that dynamics at the 360degree->0degree transition can still be identified
% INPUTS:
%   fileFilter: Search string for files to analyse, eg '*.tif'
%       Input image stacks should be cropped movies of septum/ circumferentially localized protein dynamics.
%       Cell should not overlap image edges.
%   pixSz: Camera pixel size in nanometres
% OUTPUTS:
%   Files, in subdirectory analysed:
%       <filename>_kymo.tif: Kymograph
%       <filename>_kymoWrap.tif: Kymograph, repeated nKymoWrap times
%       <filename>_fitData.mat: Parameters defining the circle used for kymograph analysis.
%       <filename>_diamInfo.txt: Ring diameter, each frame if constricting, otherwise only defined for first frame (same for the rest)
%
%
% OPTIONAL INPUTS 
%   nKymoWrap (default:2): Number of times to plot the kymograph side-by-side in the kymoWrap file.
% 
%TODO SAVE THE DIAMETER
if ~exist(nKymoWrap)
    nKymoWrap=2;
end

f = dir(fileFilter);
nF= numel(f);

hF= figure;
for ii = 1:nF
    file = f(ii).name;
    folder=f(ii).folder;
    fname = fullfile(folder,file)
    try 
        manualCircleAndKymo(fname,pixSz,nKymoWrap);
    catch ME
        display(getReport(ME));
    end
end
close(hF);
%-----------------------------------------------
function manualCircleAndKymo(fname,pixSz,nKymoWrap)


[pathstr,name,ext]=fileparts(fname);
if isempty(pathstr)
    pathstr='.';%matlab looks on all pths not just current path if dont do this
end

savepath = fullfile(pathstr,'analysed');

if ~exist(savepath,'dir')
    mkdir(savepath);
end
% calculate the kymograph
ringStack=imreadstack(fname);
[ kymo,circleData] = doManualKymo(ringStack,pixSz);


save([savepath,filesep,name,'_fitData.mat'],'circleData');
%make 2pi versions of everything for convenience
kymo_wrap = repmat(kymo,[1,nKymoWrap]);

%save everything as floats if that's what's returned
tiffwrite([savepath,filesep,name,'_kymo.tif'],kymo);
tiffwrite([savepath,filesep,name,'_kymoWrap.tif'],kymo_wrap);
