function [velocityMap, opt] = exportVelMapToVideo(velocityMap,position_t,position_x,position_y,opt)
    
if ~isfield(opt, 'exportimages') || isempty(opt.exportimages); opt.exportimages = 'n'; end
if ~isfield(opt, 'imagesformat') || isempty(opt.imagesformat); opt.imagesformat = 'png'; end
if ~isfield(opt, 'movieformat') || isempty(opt.movieformat); opt.movieformat = 'mp4'; end
if ~isfield(opt, 'OutputEvery'); opt.OutputEvery=1;end
if ~isfield(opt, 'maxV'); opt.maxV = Inf;end %((opt.ROIsize*opt.pixelSize/3)/opt.timeFrame)*60;end
if ~isfield(opt, 'timerDisplay'); opt.timerDisplay = 'UR';end
if ~isfield(opt, 'bgImage'); opt.bkgImage = 'TOI mean'; end
if strcmp(opt.bgImage, 'Original') && ~isempty(opt.filePath{1}); series1 = readFileToStack(opt.filePath{1}); end
if strcmp(opt.bgImage, 'Other') && ~isempty(opt.bg_filePath{1}); bg_series = readFileToStack(opt.bg_filePath{1}); end


alltimes=unique(position_t);

correctFactor=1;%(1/0.1)*(0.06/2);
% XLim2=[400 1100];

if strcmp(opt.movieformat,'avi')
    writerObj = VideoWriter([opt.path 'Movie_' opt.outputName '.avi'],'Uncompressed AVI');
    writerObj.FrameRate = 3;
elseif strcmp(opt.movieformat,'mp4')
    writerObj = VideoWriter([opt.path 'Movie_' opt.outputName '.mp4'],'MPEG-4');
    writerObj.FrameRate = 3;
elseif strcmp(opt.movieformat,'jpeg')
    writerObj = VideoWriter([opt.path 'Movie_' opt.outputName '.avi'],'Motion JPEG AVI');
    writerObj.FrameRate = 3;
else 
     error('Please select a valid movie format (''avi'' ,''mp4'' or ''jpeg'')');
end

open(writerObj);
%fig=figure;
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
    
    %velocityMap{i}.magnitudesPerMin = sqrt(vx(velocityMap{i}.goodVectors).^2+vy(velocityMap{i}.goodVectors).^2)*60;
%     velocityMap{i}.magnitudesPerMin = sqrt(vx.^2+vy.^2)*60; % Conversion
%     done in convertVelocitiesFromPxPy.m
    velocityMap{i}.magnitudesPerMin = velocityMap{i}.VelMap;
    
    velPerMin = velocityMap{i}.magnitudesPerMin;
    magnitudesPerMin(i,:) = velPerMin(1:size(velPerMin,1)*size(velPerMin,2));
    
    % all vector that magnitudes per minutes smaller than certain thesholdV
    badVectors = find(velocityMap{i}.magnitudesPerMin>opt.maxV);
    toigoodvect{i} = setdiff(indgoodvect,badVectors);
    allgood = union(allgood,toigoodvect{i});
    overallgoodvectors = intersect(overallgoodvectors,toigoodvect{i});
    
    velocityMap{i}.minVel = nanmin(velocityMap{i}.magnitudesPerMin(toigoodvect{i}));
    velocityMap{i}.maxVel = nanmax(velocityMap{i}.magnitudesPerMin(toigoodvect{i}));
    allminvel(i)=velocityMap{i}.minVel;
    allmaxvel(i)=velocityMap{i}.maxVel;
  
end

minGlobalVel = nanmin(nanmin(magnitudesPerMin(:,overallgoodvectors)));
maxGlobalVel = nanmax(nanmax(magnitudesPerMin(:,overallgoodvectors)));
%           
% colormapsize = length(overallgoodvectors);%
% fullcolormap = colormap(jet(colormapsize));
% fullcolormap(end,:) = [1 1 1];
% fullcolormap(1,:) = [0 0 0];
% colormap(fullcolormap);
% cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);

vectorBorder = opt.ROIsize/2;
% Only display region of image with vectors in it
% Uses same limits for all images displayed
XLim = zeros(1,2);
YLim = zeros(1,2);
         try
            XLim(1,:) = [min(posy(:)) max(posy(:))];
         catch % if there are NO good vectors in a particular map
         end
         try
            YLim(1,:) = [min(posx(:)) max(posx(:))];
         catch % if there are NO good vectors in a particular map
         end
       

for i = 1:length(velocityMap)
             
    %colormapsize = length(overallgoodvectors);%
    clear colormapsize fullcolormap
    colormapsize = length(toigoodvect{i});%
    fullcolormap = colormap(jet(colormapsize));
    fullcolormap(end,:) = [1 1 1];
    fullcolormap(1,:) = [0 0 0];
    colormap(fullcolormap);
    cMapSpacingV = (maxGlobalVel-minGlobalVel)/length(fullcolormap);

    %    fig=figure;
    %    hold on
    lowerColormapLimit = max(floor((velocityMap{i}.minVel-minGlobalVel)/cMapSpacingV),1);  

    set(gcf,'Units','normalized','Position',[.1 .1 .8 .8])
    %set(gcf,'Position',[1 1 800   800])
    subplot(10,1,1:9)
    if strcmpi(opt.bgImage, 'Original')
        imagesc(series1(:, :, ceil(position_t(i) - 0.5)));
    elseif strcmpi(opt.bgImage, 'Other')
        imagesc(bg_series(:, :, floor(position_t(i) - 0.5)));
    else %if strcmpi(opt.bgImage, 'TOI mean')
        % Code it to retrieve it from series 1 instead of saving it inside
        % of velocityMap?
        imagesc(velocityMap{i}.data_TOImean)
    end
    % imagesc(velocityMap{i}.data_myosin)
    % imagesc(velocityMap{i}.data_actin)
    
    
    if strcmp(opt.timerDisplay,'UR')
        text(0.85,0.95,[num2str(alltimes(i)*opt.timeFrame,'%1.1f') '    s'],...
        'FontSize',16,...
        'BackgroundColor',[0 0 0],...
        'Units','normalized',...
        'Color',[1 1 1]);
    elseif strcmp(opt.timerDisplay,'UL')
        text(0.05,0.95,[num2str(alltimes(i)*opt.timeFrame,'%1.1f') '    s'],...
        'FontSize',16,...
        'BackgroundColor',[0 0 0],...
        'Units','normalized',...
        'Color',[1 1 1]);
        %[0.152343750000001 0.27586909184066 0.116718750000016 0.0561563981042676],...
    elseif strcmp(opt.timerDisplay,'LL')
        text(0.05,0.05,[num2str(alltimes(i)*opt.timeFrame,'%1.1f') '    s'],...
        'FontSize',16,...
        'BackgroundColor',[0 0 0],...
        'Units','normalized',...
        'Color',[1 1 1]);
    elseif strcmp(opt.timerDisplay,'LR')
        text(0.85,0.05,[num2str(alltimes(i)*opt.timeFrame,'%1.1f') '    s'],...
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

    if ~isnan(velocityMap{i}.maxVel) % if there are no good vectors, don't plot anything!
      if cMapSpacingV == 0
         colormap(fullcolormap( lowerColormapLimit:lowerColormapLimit ,:))
      else
          colormap(fullcolormap( lowerColormapLimit:colormapsize,:))
         %colormap(fullcolormap( lowerColormapLimit:floor((velocityMap{i}.maxVel-minGlobalVel)/cMapSpacingV ) ,:))
      end
          velocityMap{i}.autoScale = quiverc(posyall(toigoodvect{i}),posxall(toigoodvect{i}),-vxall(i,toigoodvect{i}),-vyall(i,toigoodvect{i}),1);
        % velocityMap{i}.autoScale = quiverc(posyall(toigoodvect{i}),posxall(toigoodvect{i}),-vxall(i,toigoodvect{i}),-vyall(i,toigoodvect{i}),1);
         set(gca,'XTick',[],'YTick',[],'XLim',XLim,'YLim',YLim);
         fontWeight = 'normal';
         fontName = 'Arial';
         colormap(gray(colormapsize))
         hold off;

         subplot(10,1,10:10);
         subimage(ind2rgb(repmat(1:length(fullcolormap),length(fullcolormap),1),fullcolormap));
         rangeVel=minGlobalVel+(maxGlobalVel-minGlobalVel)/3:(maxGlobalVel-minGlobalVel)/3:maxGlobalVel-(maxGlobalVel-minGlobalVel)/3;
         set(gca,'XTick',[length(fullcolormap)/3 2*length(fullcolormap)/3],'XTickLabel',{num2str(rangeVel(1),'%1.2f');num2str(rangeVel(2),'%1.2f')},'YTick',[],'FontWeight',fontWeight,'FontName',fontName,'fontSize',18);

         %set(gca,'XTick',[],'YTick',[]);
         daspect([1 4.5 1])
         text(min(get(gca,'XLim')),mean(get(gca,'YLim')),{[num2str(minGlobalVel,'%1.2f') ' '];'\mum/min '},'HorizontalAlignment','right','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
         text(max(get(gca,'XLim')),mean(get(gca,'YLim')),{[' ' num2str(maxGlobalVel,'%1.2f')];' \mum/min'},'HorizontalAlignment','left','FontWeight',fontWeight,'FontName',fontName,'fontSize',20)
    % zdfsdfds
    end
    set(gcf,'Color',[1 1 1])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);

    if strcmp(opt.exportimages,'y') && mod(i,opt.OutputEvery)==0
        if strcmp(opt.imagesformat,'png')
            ext='png';
            export_fig([opt.path 'VelocityMapTime_' opt.outputName '_' num2str(round(position_t(i)*opt.timeFrame)) 's.' ext],gcf)
        elseif strcmp(opt.imagesformat,'pdf')
            ext='pdf';
            export_fig([opt.path 'VelocityMapTime_' opt.outputName '_' num2str(round(position_t(i)*opt.timeFrame)) 's.' ext],gcf)
        else
            error('Please select a valid image format (''pdf'' or ''png'')');     
        end
    end
    close              
end

close(writerObj);