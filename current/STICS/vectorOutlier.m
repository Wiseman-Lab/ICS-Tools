function [goodVectors] = vectorOutlier(gridIndexj,gridIndexi,v,subRegionStep,thresholdInput)

%Assume all vectors are good to start

goodVectors = ones(size(gridIndexj));

for i=1:length(gridIndexj)
    neighbourList = nearestNeighbour(gridIndexj(i),gridIndexi(i),gridIndexj,gridIndexi,subRegionStep,8);
    localMedian = nanmedian(v(neighbourList));
    neighbours24 = nearestNeighbour(gridIndexj(i),gridIndexi(i),gridIndexj,gridIndexi,subRegionStep,24);
    threshold = nanstd(v(neighbours24))*thresholdInput;
    if abs(localMedian - v(i)) > threshold
        goodVectors(i) = 0;
    end
end
