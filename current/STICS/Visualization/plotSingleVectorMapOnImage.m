function [velocityMap, opt] = plotSingleVectorMapOnImage(ax,velocityMap,opt, m, k)

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
            series1 = readFileToStack(opt.filePath{1,1});
            series2 = readFileToStack(opt.filePath{1,2});
        catch
            disp('Could not find the requested file, please select it manually. ');
            [fileID, filePath] = uigetfile("*.tif");
            series1 = readFileToStack([filePath, fileID]);

            [fileID, filePath] = uigetfile("*.tif");
            series2 = readFileToStack([filePath, fileID]);            
        end
    end
    if strcmp(opt.bgImage, 'Other') && ~isempty(opt.bg_filePath{1}); bg_series = readFileToStack(opt.bg_filePath{1}); end
    if ~isfield(opt, 'VelScaleRange'); opt.VelScaleRange = 'Local'; end
    
%     if ~isfield(velocityMap{m,1}, 'VelMap'); [velocityMap, opt] = convertVelocitiesFromPxPy(velocityMap, opt.position_x, opt.position_y, opt.position_t, opt); end

    alltimes=unique(opt.position_t);

    correctFactor=1;%(1/0.1)*(0.06/2);
    % XLim2=[400 1100];

%     figure;

        posx = squeeze(opt.position_x);
        posxall = posx(1:size(posx,1)*size(posx,2));
        posy = squeeze(opt.position_y);
        posyall = posy(1:size(posy,1)*size(posy,2));
    
    for i = 1:length(velocityMap)
        velocityMap{m,i}.goodVectors(end, end) = 1; % Force last vector to be valid for auto scale

        vy = squeeze(velocityMap{m,i}.vy)*correctFactor;
        vyall(i,:) = vy(1:size(vy,1)*size(vy,2));
        vx = squeeze(velocityMap{m,i}.vx)*correctFactor;
        vxall(i,:) = vx(1:size(vx,1)*size(vx,2));
        
        velocityMap{m,i}.goodVectors = reshape(velocityMap{m,i}.goodVectors, size(posx,1), size(posy,2));
        tempMap = velocityMap{m,i}.VelMap.*velocityMap{m,i}.maskedVectors;
        AllVelMap(:,i) = tempMap(1:size(tempMap,1)*size(tempMap,2));
    end
    
    switch opt.VelScaleRange
        case 'Global'
            upperVLimit = round(prctile(AllVelMap, 99.75, "all"),2);
            
            % Sneak the all times maximum into each vector field so that the auto scale
            % quiver is adjusted appropriately to compare colors accross frames.
            vxall(:,end) = -sqrt(upperVLimit);
            vyall(:,end) = -sqrt(upperVLimit);
        case 'Local'
            
            upperVLimit = round(prctile(velocityMap{m,k}.VelMap.*velocityMap{m,k}.maskedVectors, 99.75, "all"),2);
            % OR
            % upperVLimit = round(max(velocityMap{1}.VelMap,[],"all"),2);
        case 'Custom'
            upperVLimit = opt.maxV;
            
            % Sneak the all times maximum into each vector field so that the auto scale
            % quiver is adjusted appropriately to compare colors accross frames.
            vxall(:,end) = -sqrt(upperVLimit);
            vyall(:,end) = -sqrt(upperVLimit);
    end

    indgoodvect = find(velocityMap{m,1}.maskedVectors);
    overallgoodvectors = indgoodvect;
    allgood = [];

    for i = 1:length(velocityMap)
        % all vector that magnitudes per minutes smaller than certain thesholdV
        badVectors = find(velocityMap{m,i}.VelMap>upperVLimit);


        toigoodvect{i} = setdiff(indgoodvect,badVectors);
        allgood = union(allgood,toigoodvect{i});
        overallgoodvectors = intersect(overallgoodvectors,toigoodvect{i});

        velocityMap{m,i}.minVel = nanmin(velocityMap{m,i}.VelMap(toigoodvect{i}));
        velocityMap{m,i}.maxVel = nanmax(velocityMap{m,i}.VelMap(toigoodvect{i}));
        allminvel(i)=velocityMap{m,i}.minVel;
        allmaxvel(i)=velocityMap{m,i}.maxVel;
    end

%     minGlobalVel = min(min(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
%     maxGlobalVel = max(max(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
    minGlobalVel = min(min(AllVelMap, [], "omitnan"), [], "omitnan");
    maxGlobalVel = max(max(AllVelMap, [], "omitnan"), [], "omitnan");    

    opt.GlobalVelRange(m,:) = [minGlobalVel, upperVLimit, maxGlobalVel];


    


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



        %colormapsize = length(overallgoodvectors);%
        clear colormapsize fullcolormap
    %     colormapsize = length(toigoodvect{k});%
        colormapsize = length(opt.position_x)*length(opt.position_y);

        fullcolormap = colormap(ax,jet(colormapsize));
        fullcolormap(end,:) = [1 1 1];
        fullcolormap(1,:) = [0 0 0];
        colormap(ax,fullcolormap);
%         colorbar;
%         hold on;
%         caxis([0, upperVLimit]);
%         caxis('manual');
        
        
    %     cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);
%         cMapSpacingV = (upperVLimit-minGlobalVel)/length(fullcolormap);
        cMapSpacingV = upperVLimit/length(fullcolormap);

        %    fig=figure;
%            hold on
        lowerColormapLimit = max(floor((velocityMap{m,k}.minVel-minGlobalVel)/cMapSpacingV),1);  

%% Add background image
%         set(gcf,'Units','normalized','Position',[.1 .1 .8 .8])
        %set(gcf,'Position',[1 1 800   800])
%         subplot(10,1,1:9)
        if strcmpi(opt.bgImage, 'Original')
            switch m
                case 1
                    imagesc(ax, series1(:, :, ceil(opt.position_t(k) - 0.5)));
                case 2
                    imagesc(ax, series2(:, :, ceil(opt.position_t(k) - 0.5)));
                case 3
                    mixedImage = cat(3, series1(:, :, ceil(opt.position_t(k) - 0.5)), series2(:, :, ceil(opt.position_t(k) - 0.5)));
                    imagesc(ax, mean(mixedImage,3));
                case 4
                    mixedImage = cat(3, series1(:, :, ceil(opt.position_t(k) - 0.5)), series2(:, :, ceil(opt.position_t(k) - 0.5)));
                    imagesc(ax, mean(mixedImage,3));
            end
        elseif strcmpi(opt.bgImage, 'Other')
            imagesc(ax,bg_series(:, :, floor(opt.position_t(k) - 0.5)));
        elseif strcmpi(opt.bgImage, 'Black')
            imagesc(ax,zeros(size(series1,2), size(series1,1)));
        else %if strcmpi(opt.bgImage, 'TOI mean')
            % Code it to retrieve it from series 1 instead of saving it inside
            % of velocityMap?
%             imagesc(ax,velocityMap{k}.data_TOImean);
%             imagesc(ax,mean(series1(:,:,floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
           switch m
                case 1
                    imagesc(ax, mean(series1(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
                case 2
                    imagesc(ax, mean(series2(:, :, floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift)))),3));
                case 3
                    mixedImage = cat(3, mean(series1(:, :, ceil(opt.position_t(k) - 0.5)),3), mean(series2(:, :, ceil(opt.position_t(k) - 0.5)),3));
                    imagesc(ax, mean(mixedImage,3));
                case 4
                    mixedImage = cat(3, mean(series1(:, :, ceil(opt.position_t(k) - 0.5)),3), mean(series2(:, :, ceil(opt.position_t(k) - 0.5)),3));
                    imagesc(ax, mean(mixedImage,3));
            end            
        end
        
        xt = xticks; yt = yticks;
        xticklabels(ax,xt.*opt.pixelSize);
        yticklabels(ax,yt.*opt.pixelSize);        


        if strcmp(opt.timerDisplay,'UR')
            text(ax,0.85,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'UL')
            text(ax,0.05,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
            %[0.152343750000001 0.27586909184066 0.116718750000016 0.0561563981042676],...
        elseif strcmp(opt.timerDisplay,'LL')
            text(ax,0.05,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'LR')
            text(ax,0.85,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        else 
        end

        hold on
        axis image
%         title(ax,['Velocity Map for '
%         opt.axisTitle],'fontSize',18,'FontWeight','bold','Interpreter','none');
        title(ax,'Velocity Map','fontSize',18,'FontWeight','bold','Interpreter','none');
        xlabel(ax,'\mu\itm'); ylabel(ax,'\mu\itm');
        %zdfsdfds
%         colormap(ax, 'jet')

        if ~isnan(velocityMap{m,k}.maxVel) % if there are no good vectors, don't plot anything!
%           if cMapSpacingV == 0
%              colormap(ax,fullcolormap( lowerColormapLimit:lowerColormapLimit ,:));
%           else
%              colormap(ax,fullcolormap( lowerColormapLimit:colormapsize,:));
%              %colormap(fullcolormap( lowerColormapLimit:floor((velocityMap{k}.maxVel-minGlobalVel)/cMapSpacingV ) ,:))
%           end
          
          quiverc(ax, posyall(toigoodvect{k}),posxall(toigoodvect{k}),-vxall(k,toigoodvect{k}),-vyall(k,toigoodvect{k}),1);
         % velocityMap{k}.autoScale = quiverc(posyall(toigoodvect{k}),posxall(toigoodvect{k}),-vxall(k,toigoodvect{k}),-vyall(k,toigoodvect{k}),1);
%              set(gca,'XTick',[],'YTick',[],'XLim',XLim,'YLim',YLim);
         set(ax,'XLim',XLim,'YLim',YLim);
%              fontWeight = 'normal';
%              fontName = 'Arial';
%          colorbar(ax); 
         colormap(ax, gray(colormapsize));         
         
         hold off;
%% Plot colorbar
% %              subplot(10,1,10:10);
%              subimage(ind2rgb(repmat(1:length(fullcolormap),length(fullcolormap),1),fullcolormap));
%     %          rangeVel=minGlobalVel+(nn9975Percentile-minGlobalVel)/3:(nn9975Percentile-minGlobalVel)/3:nn9975Percentile-(nn9975Percentile-minGlobalVel)/3;
%              rangeVel=minGlobalVel+(upperVLimit-minGlobalVel)/3:(upperVLimit-minGlobalVel)/3:upperVLimit-(upperVLimit-minGlobalVel)/3;
%              set(gca,'XTick',[length(fullcolormap)/3 2*length(fullcolormap)/3],'XTickLabel',{num2str(rangeVel(1),'%1.2f');num2str(rangeVel(2),'%1.2f')},'YTick',[],'FontWeight',fontWeight,'FontName',fontName,'fontSize',18);
% 
%              %set(gca,'XTick',[],'YTick',[]);
%              daspect([1 4.5 1])
%              text(min(get(gca,'XLim')),mean(get(gca,'YLim')),{[num2str(minGlobalVel,'%1.2f') ' '];'\mum/min '},'HorizontalAlignment','right','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
%              text(max(get(gca,'XLim')),mean(get(gca,'YLim')),{[' ' num2str(upperVLimit,'%1.2f')];' \mum/min'},'HorizontalAlignment','left','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
%         % zdfsdfds
        end
%         set(gcf,'Color',[1 1 1])
end
