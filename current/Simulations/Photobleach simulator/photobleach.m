function [imageseries] = photobleach(matrix_final,particles,matrix1, radius, timesize, k, frames, extra)

 area = (pi*radius^2); 
 matrix = matrix_final;
 total = sum(sum(particles(:,:)));
 imageseries = zeros(size(particles,1),size(particles,2),frames+1);
imageseries(:,:,1) = matrix_final;
particles = round(matrix1./area);
 totalparticle = zeros(size(particles,1),size(particles,2));
 for it = 1:frames
     
     %
     current = sum(sum(particles(:,:)));
     %Provides number of particles expected
     intensity = round(total*exp(-k*timesize*(it)));
     nbleach = round(total*exp(-k*timesize*(it-1)))-intensity;
     %Provides a list of all the locations
     [row, col] = find(particles(:,:) > 0);
     locations = zeros(size(row,1),2);
     locations(:,1) = row;
     locations(:,2) = col;
     %Creates a list of locations with specific iterations based on number
     %of particles
     particlelocs = [];
     for i = 1:size(locations,1)
         place = particles(locations(i,1),locations(i,2));
         for j = 1:place
             particlelocs = [particlelocs; locations(i,:)];
         end
     end
     %Creating an array with random particles to remove
     list = [1:1:size(particlelocs,1)];
     list = list(randperm(length(list)));
     plocations = zeros(nbleach, 2);
     for i = 1:nbleach
         plocations(i,1) = particlelocs(list(i),1);
         plocations(i,2) = particlelocs(list(i),2);
     end
    %Creating a series of particle images
    particleimage = zeros(size(particles,1),size(particles,2));
    for i = 1:nbleach
        particleimage(plocations(i,1),plocations(i,2)) = particleimage(plocations(i,1),plocations(i,2)) + 1;
    end
    particles = particles - particleimage;
     totalparticle = totalparticle + particleimage;
     matrixtemp    = area*totalparticle;
     matrix_conv  =  round(convolveGaussian(matrixtemp,radius*6,radius));
     matrix = matrix_final - matrix_conv;
     imageseries(:,:,it+1) = matrix;
     for x = 1:size(imageseries,1)
         for y = 1:size(imageseries,2)
             if imageseries(x,y,it+1) < 0
                 imageseries(x,y,it+1) = 0;
             end
         end
     end
 end
 
 