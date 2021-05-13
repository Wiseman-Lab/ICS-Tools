function [imageseries] = meanDisp(imageseries, psf, K)

imageseries = mean(imageseries,3);
meanI = mean2(imageseries);
stdI = std2(imageseries);

for i = 1:size(imageseries,1)
    for j = 1:size(imageseries,2)
        if imageseries(i,j) > K*stdI + meanI | imageseries(i,j) < meanI - K*stdI
            imageseries(i,j) = NaN;
        end
    end
end

