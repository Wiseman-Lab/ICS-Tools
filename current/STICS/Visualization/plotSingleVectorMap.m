function [velocityMap, opt, fig, ax1, ax2] = plotSingleVectorMap(velocityMap, opt, m, k)

    if ~isfield(opt, 'axisTitle'), opt.axisTitle = [opt.fileName{1,1},...
        '_',num2str(opt.ROIsize), 'x', num2str(opt.ROIshift), '_',...
        num2str(opt.TOIsize), 'x', num2str(opt.TOIshift), '_t', num2str(opt.tauLimit),...
        'r', num2str(opt.fitRadius), 'w', num2str(opt.maxHalfWidth ),...
        'sd', num2str(opt.threshVector), 'v', num2str(opt.maxV)]; 
    end
    if ~isfield(opt, 'outputName'), opt.outputName = opt.axisTitle; end %output file name
    if ~isfield(opt, 'exportimages') || isempty(opt.exportimages); opt.exportimages = 'n'; end
    if ~isfield(opt, 'imagesformat') || isempty(opt.imagesformat); opt.imagesformat = 'png'; end
    if ~isfield(opt, 'movieformat') || isempty(opt.movieformat); opt.movieformat = 'mp4'; end
    if ~isfield(opt, 'OutputEvery'); opt.OutputEvery=1;end
    if ~isfield(opt, 'maxV'); opt.maxV = Inf; end %((opt.ROIsize*opt.pixelSize/3)/opt.timeFrame)*60;end
    if ~isfield(opt, 'timerDisplay'); opt.timerDisplay = ''; end
    if ~isfield(opt, 'bgImage'); opt.bgImage = 'Original'; end
    if strcmp(opt.bgImage, 'Original') || strcmp(opt.bgImage, 'Black')
        try
            stack1 = readFileToStack(opt.filePath{1,1});
            stack2 = readFileToStack(opt.filePath{1,2});
        catch
            disp(['Could not find ', opt.fileName{1,1}, ', please select it manually. ']);
            [fileID, filePath] = uigetfile("*.tif");
            stack1 = readFileToStack([filePath, fileID]);
            opt.filePath{1,1} = [filePath, fileID];
            if opt.STICCS
                disp(['Could not find ', opt.fileName{1,2}, ', please select it manually. ']);
                [fileID, filePath] = uigetfile("*.tif");
                stack2 = readFileToStack([filePath, fileID]);
                opt.filePath{1,2} = [filePath, fileID];
            else
                stack2 = stack1;
                opt.filePath{1,2} = [filePath, fileID];
            end
        end
    end
    if strcmp(opt.bgImage, 'Other') && ~isempty(opt.bg_filePath{1}); bg_series = readFileToStack(opt.bg_filePath{1}); end
    
    if ~isfield(opt, 'VelScaleRange') || opt.maxV == Inf
        opt.VelScaleRange = 'Local';
    elseif opt.maxV < Inf
        opt.VelScaleRange = 'Custom';
    end
    
    if ~isfield(velocityMap{m,1}, 'VelMap') && strcmpi(opt.VelUnits,'um/sec')
        [velocityMap, opt] = convertVelocitiesFromPxPy(velocityMap, opt.position_x, opt.position_y, opt.position_t, opt);
    end

    
    %%
    alltimes=unique(opt.position_t);

    correctFactor=1;%(1/0.1)*(0.06/2);% What was this for? can we remove it?
    % XLim2=[400 1100];

%     figure;

        posx = squeeze(opt.position_x);
        posxall = posx(1:size(posx,1)*size(posx,2));
        posy = squeeze(opt.position_y);
        posyall = posy(1:size(posy,1)*size(posy,2));
    
    for i = 1:length(velocityMap)
        % Force last vector to be valid for auto scale:
        velocityMap{m,i}.maskedVectors(end) = 1; 

        vy = squeeze(velocityMap{m,i}.vy)*correctFactor;
        vyall(i,:) = vy(1:size(vy,1)*size(vy,2));
        vx = squeeze(velocityMap{m,i}.vx)*correctFactor;
        vxall(i,:) = vx(1:size(vx,1)*size(vx,2));
        
        velocityMap{m,i}.goodVectors = reshape(velocityMap{m,i}.goodVectors, size(posx,1), size(posy,2));
        tempMap = velocityMap{m,i}.VelMap.*velocityMap{m,i}.maskedVectors;
        AllVelMap(:,i) = tempMap(1:size(tempMap,1)*size(tempMap,2));
    end
    
    %%
    switch opt.VelScaleRange
        case 'Global'
            upperVLimit = round(prctile(AllVelMap, 99.75, "all"),2);
            
        case 'Local'
            
%             upperVLimit = round(prctile(velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors, 99.99, "all"),2);
            % OR
            upperVLimit = round(max(velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors,[],"all"),2);
                        
        case 'Custom'
            upperVLimit = opt.maxV;
            
    end
            % Sneak the (local or all times) maximum into each vector field so that the auto scale
            % quiver is adjusted appropriately to compare colors accross frames.
            vxall(:,end) = -sqrt(upperVLimit);
            vyall(:,end) = -sqrt(upperVLimit);
            velocityMap{m,k}.VelMap(end) = upperVLimit;

%%    
    indgoodvect = find(opt.maskCell==1);
    overallgoodvectors = indgoodvect;
    allgood = []; toigoodvect = cell(size(velocityMap));
    allminvel = zeros(size(velocityMap));
    allmaxvel = zeros(size(velocityMap));
    for i = 1:length(velocityMap)
             
        toigoodvect{m,i} = find(velocityMap{m,i}.maskedVectors == 1);
        
        % all vector that magnitudes per minutes smaller than certain thesholdV           
        badVectors = find(velocityMap{m,i}.VelMap>upperVLimit);
        toigoodvect{m,i} = setdiff(toigoodvect{m,i},badVectors);
        

        allgood = union(allgood,toigoodvect{m,i});
        overallgoodvectors = intersect(overallgoodvectors,toigoodvect{m,i});

        velocityMap{m,i}.minVel = min(velocityMap{m,i}.VelMap(toigoodvect{m,i}), [], 'all', 'omitnan');
        velocityMap{m,i}.maxVel = max(velocityMap{m,i}.VelMap(toigoodvect{m,i}), [], 'all', 'omitnan');
        allminvel(i)=velocityMap{m,i}.minVel;
        allmaxvel(i)=velocityMap{m,i}.maxVel;
    end

%     minGlobalVel = min(min(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
%     maxGlobalVel = max(max(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
    minGlobalVel = min(min(AllVelMap, [], "omitnan"), [], "omitnan");
    maxGlobalVel = max(max(AllVelMap, [], "omitnan"), [], "omitnan");    

    opt.GlobalVelRange(m,:) = [minGlobalVel, upperVLimit, maxGlobalVel];


    

%%
    %           
    % colormapsize = length(overallgoodvectors);%
    % fullcolormap = colormap(jet(colormapsize));
    % fullcolormap(end,:) = [1 1 1];
    % fullcolormap(1,:) = [0 0 0];
    % colormap(fullcolormap);
    % cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);

    % Only display region of image with vectors in it
    % Uses same limits for all images displayed
    XLim = zeros(1,2);
    YLim = zeros(1,2);
             try
                XLim(1,:) = [min(posy(:))-opt.ROIshift max(posy(:))+opt.ROIshift];
             catch % if there are NO good vectors in a particular map
             end
             try
                YLim(1,:) = [min(posx(:))-opt.ROIshift max(posx(:))+opt.ROIshift];
             catch % if there are NO good vectors in a particular map
             end


%%
fig = figure; ax1 = gca;%subplot(10,1,1:9);
        %colormapsize = length(overallgoodvectors);%
        clear colormapsize fullcolormap
    %     colormapsize = length(toigoodvect{m,k});%
        colormapsize = min(length(opt.position_x)*length(opt.position_y), 2^8);

        fullcolormap = colormap(ax1, jet(colormapsize));
        fullcolormap(end,:) = [1 1 1];
        fullcolormap(1,:) = [0 0 0];
        colormap(ax1,fullcolormap);
%         colorbar;
%         hold on;
        caxis([0, upperVLimit]);
        caxis('manual');
        
        
    %     cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);
%         cMapSpacingV = (upperVLimit-minGlobalVel)/length(fullcolormap);
        cMapSpacingV = upperVLimit/length(fullcolormap);

        %    fig=figure;
%            hold on
        lowerColormapLimit = max(floor((velocityMap{m,k}.minVel-minGlobalVel)/cMapSpacingV),1);  

%% Add background image
%         set(gcf,'Units','normalized','Position',[.1 .1 .8 .8])
        %set(gcf,'Position',[1 1 800   800])
        fig.WindowState = 'maximized';
        BWcontour = bwperim(opt.dynamicMask(:,:,k), 8);
%         imshow(BWcontour);
        axis(ax1, 'ij'); hold(ax1, 'on');
        climits = [0,1200];
        
        if strcmpi(opt.bgImage, 'Original')
            switch m
                case 1
                    bgImage = stack1(:, :, ceil(opt.position_t(k) - 0.5));
%                     imagesc(ax1, stack1(:, :, ceil(opt.position_t(k) - 0.5)));
                case 2
                    bgImage = stack1(:, :, ceil(opt.position_t(k) - 0.5));
%                     imagesc(ax1, stack2(:, :, ceil(opt.position_t(k) - 0.5)));
                case 3                    
                    mixedImage = cat(3, stack1(:, :, ceil(opt.position_t(k) - 0.5)), stack2(:, :, ceil(opt.position_t(k) - 0.5)));
                    bgImage = mean(mixedImage,3);
%                     imagesc(ax1, mean(mixedImage,3));
                case 4
                    mixedImage = cat(3, stack1(:, :, ceil(opt.position_t(k) - 0.5)), stack2(:, :, ceil(opt.position_t(k) - 0.5)));
                    bgImage = mean(mixedImage,3);
%                     imagesc(ax1, mean(mixedImage,3));
            end
        elseif strcmpi(opt.bgImage, 'Other')
%             imagesc(ax1,bg_series(:, :, floor(opt.position_t(k) - 0.5)));
            bgImage = bg_series(:, :, floor(opt.position_t(k) - 0.5));
        elseif strcmpi(opt.bgImage, 'Black')
            bgImage = zeros(size(stack1,2), size(stack1,1));
%             imagesc(ax1,zeros(size(stack1,2), size(stack1,1)));
        else %if strcmpi(opt.bgImage, 'TOI mean')
            % Code it to retrieve it from series 1 instead of saving it inside
            % of velocityMap?
%             imagesc(ax1,velocityMap{k}.data_TOImean);
%             imagesc(ax1,mean(stack1(:,:,floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
           switch m
                case 1
                    bgImage = mean(stack1(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3);
%                     imagesc(ax1, mean(stack1(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
                case 2
                    bgImage = mean(stack2(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3);
%                     imagesc(ax1, mean(stack2(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
                case 3
                    mixedImage = cat(3, mean(stack1(:, :, ceil(opt.position_t(k) - 0.5)),3), mean(stack2(:, :, ceil(opt.position_t(k) - 0.5)),3));
                    bgImage = mean(mixedImage,3);
%                     imagesc(ax1, mean(mixedImage,3));
                case 4
                    mixedImage = cat(3, mean(stack1(:, :, ceil(opt.position_t(k) - 0.5)),3), mean(stack2(:, :, ceil(opt.position_t(k) - 0.5)),3));
                    bgImage = mean(mixedImage,3);
%                     imagesc(ax1, mean(mixedImage,3));
            end            
        end
        bgImage(BWcontour) = Inf;
        imagesc(ax1, bgImage, climits);
%         colorbar;
        
xticks([]); yticks([]);
%         xt = xticks; yt = yticks;
%         xticklabels(ax1,xt.*opt.pixelSize);
%         yticklabels(ax1,yt.*opt.pixelSize); 



        if strcmp(opt.timerDisplay,'UR')
            text(ax1,0.85,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'UL')
            text(ax1,0.05,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
            %[0.152343750000001 0.27586909184066 0.116718750000016 0.0561563981042676],...
        elseif strcmp(opt.timerDisplay,'LL')
            text(ax1,0.05,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'LR')
            text(ax1,0.85,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        else 
        end

        hold on
        axis image
%         title(ax1,['Velocity Map for '
%         opt.axisTitle],'fontSize',18,'FontWeight','bold','Interpreter','none');
        title(ax1,'Velocity Magnitude Map','fontSize',12,'FontWeight','bold','Interpreter','none');
%         xlabel(ax1,'\mu\itm'); ylabel(ax1,'\mu\itm');
        %zdfsdfds
%         colormap(ax1, 'jet')

        if ~isnan(velocityMap{m,k}.maxVel) % if there are no good vectors, don't plot anything!
%           if cMapSpacingV == 0
%              colormap(ax1,fullcolormap( lowerColormapLimit:lowerColormapLimit ,:));
%           else
%              colormap(ax1,fullcolormap( lowerColormapLimit:colormapsize,:));
%              %colormap(fullcolormap( lowerColormapLimit:floor((velocityMap{k}.maxVel-minGlobalVel)/cMapSpacingV ) ,:))
%           end
          
%           quiverc(ax1, posyall(toigoodvect{m,k}),posxall(toigoodvect{m,k}),-vxall(k,toigoodvect{m,k}),-vyall(k,toigoodvect{m,k}),1);
%           quiver(ax1, posyall(toigoodvect{m,k}),posxall(toigoodvect{m,k}),-vxall(k,toigoodvect{m,k}),-vyall(k,toigoodvect{m,k}),1,'m');
         % velocityMap{k}.autoScale = quiverc(posyall(toigoodvect{m,k}),posxall(toigoodvect{m,k}),-vxall(k,toigoodvect{m,k}),-vyall(k,toigoodvect{m,k}),1);
%              set(ax1,'XTick',[],'YTick',[],'XLim',XLim,'YLim',YLim);
         set(ax1,'XLim',XLim,'YLim',YLim);
             fontWeight = 'normal';
             fontName = 'Arial';
             axis ij;
             ax1.Color = 'k';
%          colorbar(ax1); 
         colormap(ax1, gray(colormapsize))         
         
         hold(ax1, 'off');
%% Plot colorbar
ax2 = [];
%              ax2 = subplot(10,1,10:10);
%              imagesc(ax2,ind2rgb(repmat(1:length(fullcolormap),length(fullcolormap),1),fullcolormap));
%     %          rangeVel=minGlobalVel+(nn9975Percentile-minGlobalVel)/3:(nn9975Percentile-minGlobalVel)/3:nn9975Percentile-(nn9975Percentile-minGlobalVel)/3;
% %              rangeVel = minGlobalVel+(upperVLimit-minGlobalVel)/3:(upperVLimit-minGlobalVel)/3:upperVLimit-(upperVLimit-minGlobalVel)/3;
%              rangeVel = [minGlobalVel, upperVLimit];
% %              set(gca,'XTick',[length(fullcolormap)/3 2*length(fullcolormap)/3],'XTickLabel',{num2str(rangeVel(1),'%1.2f');num2str(rangeVel(2),'%1.2f')},'YTick',[],'FontWeight',fontWeight,'FontName',fontName,'fontSize',18);
%              set(ax2,'XTick',[1, length(fullcolormap)],'XTickLabel',{num2str(rangeVel(1),'%1.2f');num2str(rangeVel(2),'%1.2f')},'YTick',[],'FontWeight',fontWeight,'FontName',fontName,'fontSize',8);
%              title(ax2,'Flow Velocity','fontSize',10,'FontWeight','bold','Interpreter','none');
%              
%              %set(gca,'XTick',[],'YTick',[]);
%              daspect([1 17 1])
% %              text(min(get(gca,'XLim')),mean(get(gca,'YLim')),{[num2str(minGlobalVel,'%1.2f') ' '];'\mum/min '},'HorizontalAlignment','right','FontWeight',fontWeight,'FontName',fontName,'fontSize',12)
% %              text(max(get(gca,'XLim')),mean(get(gca,'YLim')),{[' ' num2str(upperVLimit,'%1.2f')];' \mum/min'},'HorizontalAlignment','left','FontWeight',fontWeight,'FontName',fontName,'fontSize',12)
%              text(mean(get(ax2,'XLim')),max(get(gca,'YLim')),{[];'\mu\itm/min'},'HorizontalAlignment','center','FontWeight',fontWeight,'FontName',fontName,'fontSize',8);
%         % zdfsdfds
        end
%         set(gcf,'Color',[1 1 1])
end
