function batchAnalyseRing(fileFilter,kymoPixSz,kymoWidth,kymoBgFramespan,bleachPlotOn,useNewFitter,psfFWHM,useTLFitter,frameStep,doBGFG_separation)
% function batchAnalyseRing(fileFilter,kymoPixSz,kymoWidth,kymoBgFramespan,bleachPlotOn)
if ~exist('bleachPlotOn','var')
    bleachPlotOn = true;
end
if ~exist('useNewFitter','var')
    useNewFitter=false;
    psfFWHM = [];
end
if ~exist('useTLFitter','var')
    useTLFitter=false;
    frameStep=[];
end
if ~exist('doBGFG_separation','var')
    doBGFG_separation=false;
end

f = dir(fileFilter);
nF= numel(f);
for ii = 1:nF
    fname = f(ii).name;
    fname
    try 
        ringCorrectAndFit(fname,kymoPixSz,kymoWidth,kymoBgFramespan,bleachPlotOn,useNewFitter,psfFWHM,useTLFitter,frameStep,doBGFG_separation);
    catch ME
        display(getReport(ME));
    end
end

