function [outputArg1,outputArg2] = var_map (kernel_xy, black_threshold)
%VAR_MAP Summary of this function goes here
%   Detailed explanation goes here

% kernel_size = 1089;
% kernel_xy = sqrt(kernel_size);
if mod(kernel_xy, 2) == 0
    disp("Error: Kernel size must be an odd number")
end
margin = round(kernel_xy/2) - 1;

vars = zeros(x_pix - 2*margin, y_pix - 2*margin);
for x = 1:x_pix - kernel_xy
    for y = 1:y_pix - kernel_xy
        vars(x, y) = var(maxis(x : x+kernel_xy, y : y+kernel_xy, 2), 0, 'all');
    end
end

% mean(vars, 'all')
% std(vars, 0, 'all')
% histogram(vars);

ints = imbinarize(maxis(1+margin:x_pix - margin, 1+margin:y_pix - margin,1), 5000);
% imshow(ints);

if black_threshold ~= 0
    img_bw = imbinarize(vars, black_threshold);
end    
imshow(img_bw);

img_fg = maxis(1+margin:x_pix - margin, 1+margin:y_pix - margin, 2);
map = zeros(x_pix - 2*margin, y_pix - 2*margin);
for x = 1:x_pix - kernel_xy
    for y = 1:y_pix - kernel_xy
        if ints(x,y)
            if ~img_bw(x,y)
                map(x, y) = img_fg(x,y);
            end
        end
    end
end

max_map = imagesc(map);
color_map = parula(num_frames);
color_map(1,:) = 0;
colormap(color_map);

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

