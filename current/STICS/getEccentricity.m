function [eccentricityMap] = getEccentricity(velocityMap)
%getEccentricity Calculates the eccentricity of the elliptical (rotated)
%Gaussian function for the first correlation functions (tau = 0) of each
%vector.
% Purpose: The first image correlation function (tau = 0) contains
% information on the spatial distribution of features in the image. The
% first image correlation function of a single point emitter will have the
% same shape as the PSF. If the features are particle-shaped (dots or
% puncta), then the Gaussian fit will have very low eccentricity, whereas
% if the features are cable-like or elongated, then the eccentricity will
% be higher. This can be used to segment and distinguish regions of the
% image with different shapes. 
% 
% INPUT:
%       velocityMap Cell array containing structures for each vector field
%       frame. 'coeffGrotated' must be a field of each structure.
% OUTPUT:
%       eccentricityMap: A 3D array containing the value of the
%       eccentricity (between 0 and 1). 
% CREDITS:
%       Created by Rodrigo Migueles Ramirez in May 2021

    if ~isfield(velocityMap{1}, 'coeffGrotated')
        error('Field "coeffGrotated" not found within the velocityMap structure');
    else
        eccentricityMap = zeros(size(velocityMap{1}.coeffGrotated, 1), size(velocityMap{1}.coeffGrotated, 2), length(velocityMap));
        for k = 1:length(velocityMap)
            for i = 1:size(velocityMap{k}.coeffGrotated, 1)
                for j = 1:size(velocityMap{k}.coeffGrotated, 2)
                    if ~isempty(velocityMap{k}.coeffGrotated{i,j})
                        if velocityMap{k}.coeffGrotated{i,j}(1,2) >= velocityMap{k}.coeffGrotated{i,j}(1,3)
                            a = velocityMap{k}.coeffGrotated{i,j}(1,2);
                            b = velocityMap{k}.coeffGrotated{i,j}(1,3);

                        elseif velocityMap{k}.coeffGrotated{i,j}(1,2) < velocityMap{k}.coeffGrotated{i,j}(1,2)
                            b = velocityMap{k}.coeffGrotated{i,j}(1,2);
                            a = velocityMap{k}.coeffGrotated{i,j}(1,3);                
                        end

                        eccentricityMap(i,j,k) = real(sqrt(abs(1-(abs(b))^2/(abs(a))^2)));
                    else
                        eccentricityMap(i,j,k) = NaN;
                    end
                end    
            end
        end
    end
end


