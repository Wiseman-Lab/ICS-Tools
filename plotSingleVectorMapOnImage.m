function [velocityMap, opt] = plotSingleVectorMapOnImage(velocityMap,position_t,position_x,position_y, opt, k)

    if ~isfield(opt, 'axisTitle'), opt.axisTitle = [opt.fileName{1}, ...
        '_',num2str(opt.ROIsize), 'x', num2str(opt.ROIshift), '_', ...
        num2str(opt.TOIsize), 'x', num2str(opt.TOIshift), '_t', num2str(opt.tauLimit),  ...
        'r', num2str(opt.fitRadius), 'w', num2str(opt.maxHalfWidth ), ...
        'sd', num2str(opt.threshVector), 'v', num2str(opt.maxV)]; 
    end
    if ~isfield(opt, 'outputName'), opt.outputName = opt.axisTitle; end %output file name
    if ~isfield(opt, 'exportimages') || isempty(opt.exportimages); opt.exportimages = 'n'; end
    if ~isfield(opt, 'imagesformat') || isempty(opt.imagesformat); opt.imagesformat = 'png'; end
    if ~isfield(opt, 'movieformat') || isempty(opt.movieformat); opt.movieformat = 'mp4'; end
    if ~isfield(opt, 'OutputEvery'); opt.OutputEvery=1;end
    if ~isfield(opt, 'maxV'); opt.maxV = Inf; end %((opt.ROIsize*opt.pixelSize/3)/opt.timeFrame)*60;end
    if ~isfield(opt, 'timerDisplay'); opt.timerDisplay = 'LR'; end
    if ~isfield(opt, 'bgImage'); opt.bkgImage = 'TOI mean'; end
    if strcmp(opt.bgImage, 'Original') && ~isempty(opt.filePath{1})
        try
            series1 = readFileToStack(opt.filePath{1});
        catch
            disp(['Could not find file ', opt.fileName{1}]);
            [fileID, filePath] = uigetfile("*.tif");
            series1 = readFileToStack([filePath, fileID]);
        end
    end
    if strcmp(opt.bgImage, 'Other') && ~isempty(opt.bg_filePath{1}); bg_series = readFileToStack(opt.bg_filePath{1}); end


    alltimes=unique(position_t);

    correctFactor=1;%(1/0.1)*(0.06/2);
    % XLim2=[400 1100];

    figure;
        indgoodvect = find(velocityMap{1}.goodVectors);
        overallgoodvectors = indgoodvect;
        allgood = [];
    for i = 1:length(velocityMap)
        posx = squeeze(position_x);
        posxall = posx(1:size(posx,1)*size(posx,2));
        posy = squeeze(position_y);
        posyall = posy(1:size(posy,1)*size(posy,2));

        vy = squeeze(velocityMap{i}.vy)*correctFactor;
        vyall(i,:) = vy(1:size(vy,1)*size(vy,2));
        vx = squeeze(velocityMap{i}.vx)*correctFactor;
        vxall(i,:) = vx(1:size(vx,1)*size(vx,2));

    %     %velocityMap{i}.magnitudesPerMin = sqrt(vx(velocityMap{i}.goodVectors).^2+vy(velocityMap{i}.goodVectors).^2)*60;
    % %     velocityMap{i}.magnitudesPerMin = sqrt(vx.^2+vy.^2)*60; % Conversion
    % %     done in convertVelocitiesFromPxPy.m
    %     velocityMap{i}.magnitudesPerMin = velocityMap{i}.VelMap;
    %     
        tempMap = velocityMap{i}.VelMap;
        AllVelMap(:,i) = tempMap(1:size(tempMap,1)*size(tempMap,2));



    end
    nn9975Percentile = round(prctile(AllVelMap, 99.75, "all"),2);



    for i = 1:length(velocityMap)
        % all vector that magnitudes per minutes smaller than certain thesholdV
    %     badVectors = find(velocityMap{i}.VelMap>opt.thresholdV);
    %     badVectors = find(velocityMap{i}.VelMap>opt.maxV);
        badVectors = find(velocityMap{i}.VelMap>nn9975Percentile);


        toigoodvect{i} = setdiff(indgoodvect,badVectors);
        allgood = union(allgood,toigoodvect{i});
        overallgoodvectors = intersect(overallgoodvectors,toigoodvect{i});

        velocityMap{i}.minVel = nanmin(velocityMap{i}.VelMap(toigoodvect{i}));
        velocityMap{i}.maxVel = nanmax(velocityMap{i}.VelMap(toigoodvect{i}));
        allminvel(i)=velocityMap{i}.minVel;
        allmaxvel(i)=velocityMap{i}.maxVel;
    end

    minGlobalVel = min(min(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
    maxGlobalVel = max(max(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");

    opt.GlobalVelRange = [minGlobalVel, nn9975Percentile, maxGlobalVel];

    % Sneak the all times maximum into each vector field so that the auto scale
    % quiver is adjusted appropriately to compare colors accross frames.
    vxall(:,end) = -sqrt(nn9975Percentile);
    vyall(:,end) = -sqrt(nn9975Percentile);


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
        colormapsize = length(position_x)*length(position_y);

        fullcolormap = colormap(jet(colormapsize));
        fullcolormap(end,:) = [1 1 1];
        fullcolormap(1,:) = [0 0 0];
        colormap(fullcolormap);
    %     cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);
        cMapSpacingV = (nn9975Percentile-minGlobalVel)/length(fullcolormap);

        %    fig=figure;
        %    hold on
        lowerColormapLimit = max(floor((velocityMap{k}.minVel-minGlobalVel)/cMapSpacingV),1);  

        set(gcf,'Units','normalized','Position',[.1 .1 .8 .8])
        %set(gcf,'Position',[1 1 800   800])
        subplot(10,1,1:9)
        if strcmpi(opt.bgImage, 'Original')
            imagesc(series1(:, :, ceil(position_t(k) - 0.5)));
        elseif strcmpi(opt.bgImage, 'Other')
            imagesc(bg_series(:, :, floor(position_t(k) - 0.5)));
        else %if strcmpi(opt.bgImage, 'TOI mean')
            % Code it to retrieve it from series 1 instead of saving it inside
            % of velocityMap?
            imagesc(velocityMap{k}.data_TOImean)
        end


        if strcmp(opt.timerDisplay,'UR')
            text(0.85,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'UL')
            text(0.05,0.95,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
            %[0.152343750000001 0.27586909184066 0.116718750000016 0.0561563981042676],...
        elseif strcmp(opt.timerDisplay,'LL')
            text(0.05,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        elseif strcmp(opt.timerDisplay,'LR')
            text(0.85,0.05,[num2str((alltimes(k)-1)*opt.timeFrame,'%1.1f') '    s'],...
            'FontSize',16,...
            'BackgroundColor',[0 0 0],...
            'Units','normalized',...
            'Color',[1 1 1]);
        else 
        end

        hold on
        axis image
        title(['Velocity Map for ' opt.axisTitle],'fontSize',18,'FontWeight','bold','Interpreter','none')
        %zdfsdfds
        colormap('jet')

        if ~isnan(velocityMap{k}.maxVel) % if there are no good vectors, don't plot anything!
          if cMapSpacingV == 0
             colormap(fullcolormap( lowerColormapLimit:lowerColormapLimit ,:))
          else
             colormap(fullcolormap( lowerColormapLimit:colormapsize,:))
             %colormap(fullcolormap( lowerColormapLimit:floor((velocityMap{k}.maxVel-minGlobalVel)/cMapSpacingV ) ,:))
          end
             velocityMap{k}.autoScale = quiverc(posyall(toigoodvect{k}),posxall(toigoodvect{k}),-vxall(k,toigoodvect{k}),-vyall(k,toigoodvect{k}),1);
             % velocityMap{k}.autoScale = quiverc(posyall(toigoodvect{k}),posxall(toigoodvect{k}),-vxall(k,toigoodvect{k}),-vyall(k,toigoodvect{k}),1);
             set(gca,'XTick',[],'YTick',[],'XLim',XLim,'YLim',YLim);
             fontWeight = 'normal';
             fontName = 'Arial';
             colormap(gray(colormapsize))
             hold off;

             subplot(10,1,10:10);
             subimage(ind2rgb(repmat(1:length(fullcolormap),length(fullcolormap),1),fullcolormap));
    %          rangeVel=minGlobalVel+(nn9975Percentile-minGlobalVel)/3:(nn9975Percentile-minGlobalVel)/3:nn9975Percentile-(nn9975Percentile-minGlobalVel)/3;
             rangeVel=minGlobalVel+(nn9975Percentile-minGlobalVel)/3:(nn9975Percentile-minGlobalVel)/3:nn9975Percentile-(nn9975Percentile-minGlobalVel)/3;
             set(gca,'XTick',[length(fullcolormap)/3 2*length(fullcolormap)/3],'XTickLabel',{num2str(rangeVel(1),'%1.2f');num2str(rangeVel(2),'%1.2f')},'YTick',[],'FontWeight',fontWeight,'FontName',fontName,'fontSize',18);

             %set(gca,'XTick',[],'YTick',[]);
             daspect([1 4.5 1])
             text(min(get(gca,'XLim')),mean(get(gca,'YLim')),{[num2str(minGlobalVel,'%1.2f') ' '];'\mum/min '},'HorizontalAlignment','right','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
             text(max(get(gca,'XLim')),mean(get(gca,'YLim')),{[' ' num2str(nn9975Percentile,'%1.2f')];' \mum/min'},'HorizontalAlignment','left','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
        % zdfsdfds
        end
        set(gcf,'Color',[1 1 1])
end
