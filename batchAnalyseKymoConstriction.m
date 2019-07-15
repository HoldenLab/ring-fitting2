function batchAnalyseKymoConstriction(filefilter,varargin)
%function batchAnalyseKymoConstriction(filefilter)


f = dir(filefilter);
nF= numel(f)

for ii = 1:nF
    file = f(ii).name;
    folder=f(ii).folder;
    fname = fullfile(folder,file)
    fnameList{ii} = fname;
    try 
         [dt(ii),fittable,h1,diam0(ii)] = analyseKymoConstriction(fname,varargin{:});
         saveas(h1,[fname(1:end-4),'_wplot.fig']);
         saveas(h1,[fname(1:end-4),'_wplot.png']);
         writetable(fittable,[fname(1:end-4),'_wdata.csv'])
    catch ME
        display(getReport(ME));
    end
end

t=table(fnameList',dt',diam0','VariableNames',{'Filename','ConstrictionTime','InitialDiameter'});
savename= fullfile(folder,'constrictionTimeList.csv');
writetable(t,savename);

