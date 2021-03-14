function imageMatrix=createImage(sizeImage, numParticles)

% creates a matrix with numParticles number of particles 

xCoor=ceil(sizeImage*rand(numParticles, 1));
yCoor=ceil(sizeImage*rand(numParticles, 1));

imageMatrix=full(sparse(xCoor, yCoor, 1, sizeImage, sizeImage));
