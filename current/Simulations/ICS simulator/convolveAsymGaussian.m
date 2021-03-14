function convolvedMatrix=convolveAsymGaussian(matrix, filterSize, radius)




standardDev=radius/2;
v = fspecial( 'gaussian', [filterSize 1], standardDev(1)); % vertical filter
h = fspecial( 'gaussian', [1 filterSize], standardDev(2) ); % horizontal

filter=v * h;
convolvedMatrix=imfilter(double(matrix),filter,'circular','conv');