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

