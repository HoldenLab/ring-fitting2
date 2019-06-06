function nms = nmsfilter(im,nms_con)


%3.non-maximum suppression
[nR nC] = size(im);
nms = zeros(nR, nC);
for ii = 2:nR-1 %could make this better by padding
    for jj=2:nC-1
        if strcmp(nms_con,'horz')%this makes the most sense for kymographs
            if (im(ii,jj) == max([im(ii,jj), im(ii,jj+1), im(ii,jj-1)])) 
                nms(ii,jj) = im(ii,jj);
            end
        elseif strcmp(nms_con,'4-con')
            if (im(ii,jj) == max([im(ii,jj), im(ii+1,jj), im(ii-1,jj)])) || ...
                (im(ii,jj) == max([im(ii,jj), im(ii,jj+1), im(ii,jj-1)])) 
                    nms(ii,jj) = im(ii,jj);
            end
        elseif strcmp(nms_con,'8-con')
            if (im(ii,jj) == max([im(ii,jj), im(ii+1,jj+1), im(ii-1,jj-1)])) || ...
                (im(ii,jj) == max([im(ii,jj), im(ii-1,jj+1), im(ii+1,jj-1)]))
                    nms(ii,jj) = im(ii,jj);
            end
        else
            error('nmsfilter:unrecognizedNmsArg','Unregonized NMS connectivity argument');
        end
    end
end

