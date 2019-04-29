function [fg bg] = kymoMedFilter(kymoIm, frameSpan)
%returns double

nCol = size(kymoIm,2);
bg = zeros(size(kymoIm));
for ii = 1:nCol
    bg(:,ii) = medfilt1(kymoIm(:,ii),frameSpan,'includenan','truncate');
end
fg = kymoIm-bg;

