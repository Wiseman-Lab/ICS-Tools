function convolvedMatrix=convolveDisk(matrix, radius)

%uses built in Matlab filer function to convolve particle postions (matrix) with disk like PSF
%h = fspecial('disk', radius) returns a circular averaging filter (pillbox) 
%within the square matrix of side 2*radius+1. The default radius is 5.

filter=fspecial('disk', radius);
convolvedMatrix=imfilter(double(matrix),filter,'circular','conv');