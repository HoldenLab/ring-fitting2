function trackResults = appendDiameter(fname)
%function trackResults = appendDiameter(fname)
% read the diameters from each mat file and append them to the analysed
% track speeds
% saves the data as a new csv file
% optionally outputs the data as a table to matlab
trackResults = readtable(fname);

%diameter list
f= dir('*_fitData.mat');
diamData = table();
tmpTable=table();
for ii = 1:numel(f)
    tmpTable.FileName =  {f(ii).name};
    data = importdata(f(ii).name);
    tmpTable.DiameterNm = data.kymoInfo(1,2)*2;
    diamData=[diamData;tmpTable];
end

DiameterNm = zeros(size(trackResults,1),1);
trackResults.DiameterNm = DiameterNm;


nFile = size(diamData,1);
for ii = 1:nFile
    fCur = diamData.FileName{ii};
    fStub = strrep(fCur,'_fitData.mat','');
    isCurFile = find(startsWith(trackResults.Image_ROI_Name,[fStub,'_']));
    trackResults.DiameterNm(isCurFile)= diamData.DiameterNm(ii);
end

writetable(trackResults,[fname(1:end-4),'_diam.csv'])
