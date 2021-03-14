function [population] = simul8trMovement3D(population,timesize,pixelsize,sizeX,sizeY,sizeZ);


    population.xCoor = population.xCoor+(population.flowX*timesize)/pixelsize;
    population.yCoor = population.yCoor+(population.flowY*timesize)/pixelsize;0
    population.zCoor = population.zCoor+(population.flowZ*timesize)/pixelsize;


% Fix BC problems here
population.xCoor = mod(population.xCoor,sizeX);
population.yCoor = mod(population.yCoor,sizeY);
population.xCoorDisplay = round(population.xCoor);
population.yCoorDisplay = round(population.yCoor);
xZeros = find(population.xCoorDisplay==0);
yZeros = find(population.yCoorDisplay==0);
population.xCoorDisplay(xZeros) = sizeX;
population.yCoorDisplay(yZeros) = sizeY;
if sizeZ ~= 0
    population.zCoor = mod(population.zCoor,sizeZ);
end