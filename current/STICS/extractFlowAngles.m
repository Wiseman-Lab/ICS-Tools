function [velocityMap, opt] = extractFlowAngles(velocityMap, opt)



for m = 1:size(velocityMap,1)
    for k = opt.kRange
        
        velocityMap{m,k}.vx = -velocityMap{m,k}.vx;
        
        TanOmega = velocityMap{m,k}.vy./velocityMap{m,k}.vx;
        Theta = atan(TanOmega);
        
        xPyP = velocityMap{m,k}.vx>0 & velocityMap{m,k}.vy>0;
        Omega = xPyP.*Theta;
        
        xPyN = velocityMap{m,k}.vx>0 & velocityMap{m,k}.vy<0;
        Omega = Omega + xPyN.*(2*pi + Theta);
        
        xNyP = velocityMap{m,k}.vx<0 & velocityMap{m,k}.vy>0;
        Omega = Omega + xNyP.*(pi + Theta);
        
        xNyN = velocityMap{m,k}.vx<0 & velocityMap{m,k}.vy<0;
        Omega = Omega + xNyN.*(pi + Theta);
        
        
        velocityMap{m,k}.vx = -velocityMap{m,k}.vx;
        
        if ~isfield(opt, 'AngleUnits')
            opt.AngleUnits = 'Radians';
        elseif strcmpi(opt.AngleUnits, 'Degrees')
            Omega = rad2deg(Omega);
        end            
        
        velocityMap{m,k}.FlowAngles = Omega;
        
%         Theta(:,:,k) = atand(velocityMap{k}.vy./-velocityMap{k}.vx);        
%         for i = 1:size(position_x, 1)
%             for j = 1:size(position_y, 2)
%                 if velocityMap{k}.vx(i,j)>0 && velocityMap{k}.vy(i,j)>0
%                     Omega(i,j,k) = Theta(i,j,k);
%                 elseif velocityMap{k}.vx(i,j)<0 && velocityMap{k}.vy(i,j)>0
%                     Omega(i,j,k) = 180 - (-Theta(i,j,k));
%                 elseif velocityMap{k}.vx(i,j)>0 && velocityMap{k}.vy(i,j)<0
%                     Omega(i,j,k) = 360 - (-Theta(i,j,k));
%                 elseif velocityMap{k}.vx(i,j)<0 && velocityMap{k}.vy(i,j)<0
%                     Omega(i,j,k) = 180 + Theta(i,j,k);
%                 end
%             end
%         end
    end

% figure; imagesc(velocityMap{m,k}.FlowAngles); axis square; colormap(hsv);
% c = colorbar; c.Label.String = ['Angle in ', opt.AngleUnits]; caxis([0 2*pi]);
% hold on;
% quiver(repmat((1:size(opt.vectorPositions,1)), size(opt.vectorPositions,1), 1),...
%     repmat((1:size(opt.vectorPositions,1))', 1 ,size(opt.vectorPositions,1)),...
%     -velocityMap{m,k}.vx, -velocityMap{m,k}.vy,1, 'k', 'LineWidth', 1);

end




%     
% figure; 
% % ax = subplot(1,2,1);
% % quiverc(ax, position_x, position_y, velocityMap{1,1}.vx, velocityMap{1,1}.vy,1); 
% % stack = readFileToStack(opt.filePath{1,1});
% % xlim([-size(stack,1), size(stack,1)]);
% % ylim([-size(stack,2), size(stack,2)]);
% % set(ax, "Color", "k");
% % axis square
% % 
% % ax =subplot(1,2,2);
% ax = gca;
% quiverc(ax, position_y, position_x, -velocityMap{m,k}.vx, -velocityMap{m,k,1}.vy, 1); 
% % xlim([-size(stack,1), size(stack,1)]);
% % ylim([-size(stack,2), size(stack,2)]);
% set(ax, "Color", "k");
% axis square ij
% 
% 
% 
% figure; imagesc(Omega); axis square; colormap(hsv); colorbar; caxis([0 2*pi]);
% 
% figure; polarhistogram(Omega);
% 
% figure; compass(-velocityMap{m,k}.vx, velocityMap{m,k}.vy);
% 
% 
% % max(velocityMap{1,1}.VelMap)
% % [V, I] = max(velocityMap{1,1}.VelMap,[], "all")
% % [V, I] = max(velocityMap{1,1}.VelMap,[], "all", "linear")
% % position_x(I)
% % position_y(I)
% % velocityMap{1,1}.VelMap(I)
% % set(ax,"Color", "k");
% % velocityMap{1,1}.vx(I)
% % velocityMap{1,1}.vy(I)
% % figure; ax=gca; quiverc(ax, position_x(I), position_y(I), -velocityMap{1,1}.vx(I), -velocityMap{1,1}.vy(I));
% % set(ax, "Color", "k");
% % colorbar;
% % figure; ax=gca; quiver(ax, position_x(I), position_y(I), -velocityMap{1,1}.vx(I), -velocityMap{1,1}.vy(I));
% % axis image
% % axis ij
% % figure; ax=gca; quiver(ax, position_y(I), position_x(I), -velocityMap{1,1}.vx(I), -velocityMap{1,1}.vy(I));
% % axis image ij
% % -velocityMap{1,1}.vy(I)
% % figure; ax=gca; quiver(ax, position_x(I), position_y(I), velocityMap{1,1}.vx(I), velocityMap{1,1}.vy(I)); axis image ij;
% % T = velocityMap{1,1}.vy(I)/velocityMap{1,1}.vx(I)
% % alpha = atan(T)
% % rad2deg(alpha)
% % velocityMap{1,1}.vy(I)
% % -velocityMap{1,1}.vy(I)/-velocityMap{1,1}.vx(I)
% % 
% % figure; compass(-velocityMap{1,1}.vx, velocityMap{1,1}.vy);
    
end