function [fig, ax, cb, q] = plotFlowAngles(velocityMap, opt, m, k, dark_mode, showQuiver, interpolateVectors)
    
    cmap = VectorColormap(10, 360, dark_mode);    
    
    fig = figure;  fig.WindowState = 'maximized'; ax = gca;
    
    if interpolateVectors
        vX = imresize(velocityMap{m,k}.vx, opt.ROIshift, 'bilinear');
        vY = imresize(velocityMap{m,k}.vy, opt.ROIshift, 'bilinear');
    else
        vX = velocityMap{m,k}.vx;
        vY = velocityMap{m,k}.vy;
    end
    
    if isfield(velocityMap{m,k}, 'FlowAngles')
        angles = velocityMap{m,k}.FlowAngles;
    else
        [velocityMap, opt] = extractFlowAngles(velocityMap, opt);
        angles = velocityMap{m,k}.FlowAngles;
    end
    
    angles(isnan(angles)) = 0;
    velocities = velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors;    
    
    
%     if strcmpi(opt.AngleUnits, 'Radians')
%         angles(angles <= 0.0246) = 2*pi;
%     else
%         angles(angles <= 1.4063) = 360;
%     end
%     angles(velocities <= 0.01) = 0;

    if opt.maxV == Inf
        upperVLimit = round(prctile(velocities, 99.75, "all"),2);
    else
        upperVLimit = opt.maxV;
    end
    lowerVLimit = 0.01;%round(prctile(velocityMap{m,k}.VelMap, 10, "all"),2);    
    nVelocities = velocities./upperVLimit;
    nVelocities(nVelocities < lowerVLimit/upperVLimit) = 0;
    nVelocities(nVelocities > 1) = 0;
    nVelocities(isnan(nVelocities)) = 0;
    binSize = 1/size(cmap,3);
    levels = round(nVelocities./binSize);   
    
    degrees = strcmpi(opt.AngleUnits, 'Degrees');
    if degrees
        nAngles = round(angles./360).*size(cmap,1)/10;
    else
        nAngles = round(angles./2*pi).*size(cmap,1)/10;
    end
    nAngles(velocityMap{m,k}.maskedVectors == 0) = NaN;
    nAngles(nAngles == 0) = 360;

    imFlow = zeros(size(levels,1), size(levels,2), 3);

    for ix = 1:size(imFlow,1)
        for iy = 1:size(imFlow,2)
            if levels(ix,iy) ~= 0
                imFlow(ix, iy, :) = cmap(nAngles(ix, iy),:,levels(ix, iy));
            elseif dark_mode
                imFlow(ix, iy, :) = [0, 0, 0];
            elseif ~dark_mode
                imFlow(ix, iy, :) = [1, 1, 1];
            end
        end
    end

    image(ax, imFlow); axis(ax, 'image');
    
    xlabel('Vector index (X)'); ylabel('Vector index (Y)');
    colormap(ax, cmap(:,:,end));
    cb = colorbar; cb.Label.String = ['Angle in ', opt.AngleUnits]; 
    if degrees; caxis([0, 360]); else; caxis([0 2*pi]); end
    
    hold on;
    
    nonZeroVelocities = zeros(size(nVelocities));
    nonZeroVelocities(nVelocities ~= 0) = 1;
    
    if showQuiver
        q = quiver(repmat((1:size(opt.maskedPositions,2)), size(opt.maskedPositions,1), 1),...
            repmat((1:size(opt.maskedPositions,1))', 1 ,size(opt.maskedPositions,2)),...
            -vX.*nonZeroVelocities,...
            -vY.*nonZeroVelocities,...
            1, 'LineWidth', 1);
        if dark_mode
            q.Color = 'w';
        else
            q.Color = 'k';
        end
    else
        q = [];
    end    
end

% %%
% % Script
% u = velocityMap{1,1}.vx;
% v = velocityMap{1,1}.vy;
% degrees = strcmpi(opt.AngleUnits, 'Degrees');
% % Omega = velocityMap{1,1}.FlowAngles;
% Omega = getAngles(u, v, degrees);
% cmap = VectorColormap(depth, n_colors, dark_background);
% % figure; plotAnglesImage(x, y, u, v, cmap);
% 
% figure; histogram(velocityMap{1,1}.VelMap);
% upperVLimit = round(prctile(velocityMap{m,k}.VelMap, 99.75, "all"),2);
% lowerVLimit = 0.01;%round(prctile(velocityMap{m,k}.VelMap, 10, "all"),2);
% nVelocities = velocityMap{1,1}.VelMap./upperVLimit;
% outBound = nVelocities(nVelocities > 1);
% nVelocities(nVelocities > 1) = 0;
% binSize = 1/depth;
% levels = round(nVelocities./binSize);
% levels(levels==0) = 1;
% if degrees
%     nAngles = round(Omega./360).*n_colors/10;
% else
%     nAngles = round(Omega./2*pi).*n_colors/10;
% end
% nAngles(nAngles == 0) = 360;
% 
% flowImage = zeros(size(levels,1), size(levels,2), 3);
% 
% for ix = 1:size(flowImage,1)
%     for iy = 1:size(flowImage,2)
%         if levels(ix,iy) ~= 0
%             flowImage(ix, iy, :) = cmap(nAngles(ix, iy),:,levels(ix, iy));
%         elseif dark_background
%             flowImage(ix, iy, :) = [0, 0, 0];
%         elseif ~dark_background
%             flowImage(ix, iy, :) = [1, 1, 1];
%         end
%     end
% end

