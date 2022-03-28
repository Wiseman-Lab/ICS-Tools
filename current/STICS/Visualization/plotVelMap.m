function [fig, ax, c] = plotVelMap(velocityMap, opt, m, k)
    fig = figure; fig.WindowState = 'maximized'; ax = gca;
        
    BWcontour = bwperim(opt.dynamicMask(:,:,floor(opt.position_t(k))), 8);    
    
    bwparula = colormap('parula'); 
    bwparula(1,:) = [0, 0, 0];
    bwparula(end,:) = [1, 1, 1];
    colormap(bwparula);                   
    
    
    try
        if opt.Interpolate            
            canvas = zeros(size(BWcontour));
            if isfield(velocityMap{m,k}, 'iVelMap')
                VelMapImage = velocityMap{m, k}.iVelMap.*opt.dynamicMask(:,:,floor(opt.position_t(k)));
                VelMapImage(VelMapImage>opt.maxV) = Inf;
            else
                VelMapImage = velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors;
                VelMapImage = imresize(VelMapImage, opt.ROIshift, 'bilinear');
                VelMapImage(VelMapImage>opt.maxV) = Inf;
            end
            sd = size(BWcontour)-size(VelMapImage);
            
            if sd(1)~=0 || sd(2)~=0
                canvas(ceil(sd(1)/2):floor(sd(1)/2)+size(VelMapImage,1),...
                    ceil(sd(2)/2):floor(sd(2)/2)+size(VelMapImage,2)) = VelMapImage;
            else
                canvas = VelMapImage;
            end
            canvas(BWcontour) = max(VelMapImage, [], 'all');
%             canvas(~opt.dynamicMask(:,:,opt.position_t(K))) = 0;
            imagesc(ax, canvas); 
        else
            VelMapImage = velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors;
            VelMapImage(VelMapImage>opt.maxV) = Inf;            
            imagesc(ax, VelMapImage); 
        end
    catch
        VelMapImage = velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors;
        VelMapImage(VelMapImage>opt.maxV) = Inf;         
        imagesc(ax, VelMapImage); 
    end
    
%     if scaled_axis
%         xt = xticks; yt = yticks;
%         xticklabels(xt.*opt.pixelSize);
%         yticklabels(yt.*opt.pixelSize);
%         xlabel('Distance in X (um)'); 
%         ylabel('Distance in Y (um)'); 
%     else
        xlabel('Vector index (X)'); ylabel('Vector index (Y)');
%     end
    
    c = colorbar; c.Label.String = ['Flow velocity (', opt.VelUnits, ')']; 
    axis(ax, 'ij'); axis(ax, 'image'); 
%     hold on;
%     q = quiver(repmat((1:size(opt.vectorPositions,2)), size(opt.vectorPositions,1), 1),...
%         repmat((1:size(opt.vectorPositions,1))', 1, size(opt.vectorPositions,2)),...
%         -velocityMap{m,k}.vx, -velocityMap{m,k}.vy,1, 'k', 'LineWidth', 1);
% %     q = quiver(repmat((1:size(velocityMap{m,k}.goodVectors,2)), size(velocityMap{m,k}.goodVectors,1), 1),...
% %         repmat((1:size(velocityMap{m,k}.goodVectors,1))', 1, size(velocityMap{m,k}.goodVectors,2)),...
% %         -velocityMap{m,k}.vx, -velocityMap{m,k}.vy,1, 'k', 'LineWidth', 1);
end