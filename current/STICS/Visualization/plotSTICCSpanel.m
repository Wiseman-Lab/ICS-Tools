function [velocityMap, opt] = plotSTICCSpanel(velocityMap, position_t, position_x, position_y, opt, k)

    
    for m = 1:4 % Corr11, Corr22, Corr12, Corr21
        sps = [1,4,2,3]; % Subplot indexes
        ax = subplot(2,2,sps(m));
        [velocityMap, opt] = plotSingleVectorMapOnImage(ax, velocityMap, position_t, position_x, position_y, opt, k, m);        
    end
    
    for m = 1:4
        ax = subplot(2,2,m); 
        colormap(ax, gray);        
    end
    
end