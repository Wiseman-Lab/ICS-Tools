function cmap = VectorColormap(depth, n_colors, dark_background)
%% Custom Vector Colormap
% Generates a custom colormap based on HSV, that can be used to label
% vectors based on their properties: The angle determines the hue (color)
% of the vector, while the size or magnitude of the angle determines either
% its saturation (from white (0) to full-color (1)) or its value (from
% black (0) to full-color (1)), depending on the logical variable dark_background.
% Created by Rodrigo Migueles Ramirez. McGill University. February 2022.

% n_colors = 360;
% depth = 10;
cmap = zeros(n_colors, 3, depth);
% dark_background = 0;

% figure;
for k = 1:depth %Depth
    rgbc = hsv(n_colors); % Generate HSV colormap in RGB space
    hsvc = rgb2hsv(rgbc); % convert to HSV color space
    switch dark_background
        case 0
            hsvc(:,2) = k/depth; % Modulate the Saturation
            cmap(:,:,k) = hsv2rgb(hsvc);            
        case 1
            hsvc(:,3) = k/depth; % Modulate the Value
            cmap(:,:,k) = hsv2rgb(hsvc);
    end
%     subplot(depth,1,k); imagesc(1:360); ax = gca; colormap(ax,cmap(:,:,k));
end




