fname='180327_mNG-ftsz1.tif';
fname2='180327_mNG-ftsz1_manual.tif';
pixSzNm = 65; % [nm];

copyfile(fname,fname2);
manualVerciniAnalysis(fname2,pixSzNm);
delete(fname2);
