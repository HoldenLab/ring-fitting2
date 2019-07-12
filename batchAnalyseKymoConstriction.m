function batchAnalyseKymoConstriction(filefilter)
%function batchAnalyseKymoConstriction(filefilter)


f = dir(filefilter);
nF= numel(f)

for ii = 1:nF
    file = f(ii).name;
    folder=f(ii).folder;
    fname = fullfile(folder,file)
    fnameList{ii} = fname;
    try 
         [dt(ii),fittable,h1] = analyseKymoConstriction(fname);
         saveas(h1,[fname(1:end-4),'_wplot.fig']);
         saveas(h1,[fname(1:end-4),'_wplot.png']);
         writetable(fittable,[fname(1:end-4),'_wdata.csv'])
    catch ME
        display(getReport(ME));
    end
end

t=table(fnameList',dt','VariableNames',{'Filename','ConstrictionTime'});
savename= fullfile(folder,'constrictionTimeList.csv');
writetable(t,savename);

