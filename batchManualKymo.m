function batchManualKymo(fileFilter,pixSz)
NKYMOWRAP=3;

f = dir(fileFilter);
nF= numel(f);

hF= figure;
for ii = 1:nF
    file = f(ii).name;
    folder=f(ii).folder;
    fname = fullfile(folder,file)
    try 
        manualCircleAndKymo(fname,pixSz,NKYMOWRAP);
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
