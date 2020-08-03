function nearestNeighbourList = nearestNeighbour(pointX,pointY,listX,listY,dist,numNeighbours)

% Given a point located at (pointX,pointY) finds its nearest 8 neighbours,
% assuming a distance dist given in pixels between adjacent points.  Returns
% nearestNeighbourList, which are the indices for listX and listY of the
% neighbours.
% By DK, Jan 8, 2008
% Usage:
% nearestNeighbourList =
% nearestNeighbour(pointX,pointY,listX,listY,distance,numNeighbours);

if (numNeighbours ~= 8)&&(numNeighbours ~= 24)
    error('numNeighbours must be either 8 or 24');
end

% Use a threshold for comparing floating point numbers to avoid round-off
% error
floatingThresh = eps*1e5;

nearestNeighbourList = [];

distanceList(1) = dist;
distanceList(2) = sqrt(2*dist^2);

if numNeighbours == 24
    distanceList(3) = dist*2;
    distanceList(4) = sqrt(2*(2*dist)^2);
    distanceList(5) = sqrt((2*dist)^2+dist^2);
end

for i = 1:length(distanceList);
    nearestNeighbourList = vertcat(nearestNeighbourList,find((abs(distance([pointX pointY]',[listX listY]')-distanceList(i)))<floatingThresh)');
end