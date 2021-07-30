function [velocityMap, position_x, position_y, position_t, opt] = stics_vectormapping_2(series, opt)
%this functions takes the timeseries and calculates for defined regions of
%size 'ROIsize' and 'TOIsize' stics correlation function and fits symmetric
%Gaussian to extract the flow in x and y directions...
%Developed by Elvis Pandzic 2012-2014 under supervision of 
%Prof. Paul W. Wiseman @ McGill University. 
% correspondence addressed to : e.pandzic@unsw.edu.au


% Look at opt structure and either assign default value, or use passed
% value
if nargin==1
    opt.pixelSize = 0.1; % just initialize it to something so the code below doesn't throw an error
end
if ~isfield(opt, 'pixelSize'), opt.pixelSize = 1; end %pixel size in um
if ~isfield(opt, 'timeFrame'), opt.timeFrame = 1; end %time delay (s) between subsequent frames
if ~isfield(opt, 'tauLimit'), opt.tauLimit = 130; end % what is highest time lag to compute stics of TO
%filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
if ~isfield(opt, 'filtering'), opt.filtering = 'none'; end
if ~isfield(opt, 'MoveAverage'), opt.MoveAverage = 21; end
if ~isfield(opt, 'fitRadius'), opt.fitRadius = 32; end %how many pixels (radius) are considered around the peak of corr fn when fitting 2D Gaussian
if ~isfield(opt, 'omegaThreshold'), opt.omegaThreshold = 32; end % % how big (pixels) do you allow the radius of corr fn to be...removes wide corr fns...
if ~isfield(opt, 'threshVector'), opt.threshVector = 5; end % what is the threshold delta(v) (um/s) used in discarding spurious vectors
%if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
if ~isfield(opt, 'ROIsize'), opt.ROIsize = 32; end % ROI size in pixels
if ~isfield(opt, 'ROIshift'), opt.ROIshift = 8; end %what is ROI centers shift in pixels
if ~isfield(opt, 'TOIsize'), opt.TOIsize = 5; end % what is TOI size in frames
if ~isfield(opt, 'TOIshift'), opt.TOIshift = 1; end %how much is TOI shifted
if ~isfield(opt, 'axisTitle'), opt.axisTitle = [opt.fileName{1}, ...
        num2str(opt.ROIsize), 'x', num2str(opt.ROIshift), '_', ...
        num2str(opt.TOIsize), 'x', num2str(opt.TOIshift), ...
        'Rd', num2str(opt.fitRadius), 'Om', num2str(opt.omegaThreshold ), 'Va', num2str(opt.threshVector)]; 
end %what will be used on vector map title
if ~isfield(opt, 'outputName'), opt.outputName = opt.axisTitle; end %output file name
if ~isfield(opt, 'path'), opt.path = cd; end
if ~isfield(opt, 'exportimages'), opt.exportimages ='n'; end %'y' if you want export a pdf of every vector map, 'n' if you do not
if ~isfield(opt, 'imagesformat'), opt.imagesformat ='png'; end % 'png' or 'pdf' images...can add other option (formats) in plotSingleVectorMap code
if ~isfield(opt, 'movieformat'), opt.movieformat ='mp4'; end % movie format can be avi, jpeg or mp4


% shift of ROI's amount & define the position of ROI-TOI
fracROIshift = opt.ROIsize/opt.ROIshift;
fracTOIshift = opt.TOIsize/opt.TOIshift;

along_y = floor((size(series,2)-opt.ROIsize)/opt.ROIshift)+1;
along_x = floor((size(series,1)-opt.ROIsize)/opt.ROIshift)+1;
along_t = floor((size(series,3)-opt.TOIsize)/opt.TOIshift)+1;

position_x  =   zeros(1,along_x);
position_y  =   zeros(1,along_y);
position_t  =   zeros(1,along_t);

% defining all of the ROI centers positons in x and y
for i=1:along_x
    position_x(i) = 1/2 + (i-1)/fracROIshift*opt.ROIsize + opt.ROIsize/2;
end
position_x=repmat(position_x,along_y,1);
position_x=position_x';

for j=1:along_y
    position_y(j) = 1/2 + (j-1)/fracROIshift*opt.ROIsize + opt.ROIsize/2;
end
position_y=repmat(position_y,along_x,1);

% %pre-process to determine if you want to analyse the whole FOV or select a
% %polygon to analyze...in case of polygon, it will speed up the whole
% %calculations because less ROI-TOI to consider...
% imagesc(mean(series,3)) 
% prompt = {'want to select the polygon or use whole FOV?'};
%     dlg_title = 'y or n';
%     num_lines = 1;
%     def = {'y'};
%     answer = inputdlg(prompt,dlg_title,num_lines,def);
% if strcmp(answer,'y')
%     maskCell=[];
%     while isempty(maskCell)
%         [x,y,maskCell,poly1.x,poly1.y]=roipoly;
%     end
% elseif strcmp(answer,'n')
%     maskCell=ones(size(series,1),size(series,2));
% end
% close
% opt.maskCell = maskCell;

vectorPositions = zeros(along_x,along_y);
for w=1:along_x
    for s=1:along_y
        indx = 1+(w-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((w-1)/fracROIshift)));
        indy = 1+(s-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((s-1)/fracROIshift))); 
        indx = floor(indx);
        indy = floor(indy);
        if(opt.maskCell(indx,indy)==1)
            vectorPositions(w,s)=1;
        end
    end
end
opt.vectorPositions = vectorPositions;

%immobile filter data if necessary
if strcmp(opt.filtering,'FourierWhole')
    series = immfilter(series,'F',1); 
elseif strcmp(opt.filtering,'MovingAverage')
    series = immfilter(series,opt.MoveAverage); 
elseif strcmp(opt.filtering,'butterIIR')
    %set the cutoff frequency to be the 2*sampling frequency / size of the series 
    % this way we will certainly remove the immobile component
    series = butterIIR(series, 2/(opt.timeFrame*size(series,3)),'none');
elseif strcmp(opt.filtering,'none')
else
    error('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
end

% %indexes of ROI within the mask of interest (polygon)
% [iroi, jroi] = find(vectorPositions);
velocityMap = cell(1, along_t);

parfor k=1:along_t % loop along time
    [iroi, jroi] = find(vectorPositions);
    options = opt;
    stack = series;
    % position of TOI centers in time
    position_t(k) = 1/2 + (k-1)/fracTOIshift*options.TOIsize + options.TOIsize/2;
    % define whole FOV TOI
    
    TOIFOV = stack(:,:,floor(1+(k-1)/fracTOIshift*options.TOIsize):floor(options.TOIsize*(1+((k-1)/fracTOIshift))));
    % define the average FOV TOI image that will be used in the display of the vector maps
    velocityMap{k}.data_TOImean = mean(TOIFOV, 3);
    
    
    
    tStart = tic; % 
    fprintf('analyzing TOI %i of %i ____ date and time is %s',k,along_t,datestr(now)); 
    
    velocityMap{k}.sigFitRatio = zeros(along_x, along_y);
    velocityMap{k}.coeffGtime = cell(along_x, along_y);
    velocityMap{k}.Vx = zeros(along_x, along_y, 2);
    velocityMap{k}.Vy = zeros(along_x, along_y, 2);
    velocityMap{k}.VelMap = zeros(along_x, along_y);
    
    for u=1:length(iroi) %accessing the data of only 
     %the ROI that are defined within the selected polygon
         i=iroi(u);
         j=jroi(u);
            %define region (x,y,t) on which stics will be applied
            regionanalyse = TOIFOV(floor(1+(i-1)/fracROIshift*options.ROIsize):floor(options.ROIsize*(1+((i-1)/fracROIshift))), floor(1+(j-1)/fracROIshift*options.ROIsize): floor(options.ROIsize*(1+((j-1)/fracROIshift))),:);
           
            %apply regular cross-corr stics
            [corrfn] = stics(regionanalyse, options.tauLimit);
            velocityMap{k}.sigFitRatio(i,j) = size(corrfn,3)/size(regionanalyse,3);
            %[i j k]
         
            %fit vxtau and vytau vs tau linearly to extract vx and vy
            if size(corrfn,3)>1 && size(corrfn,1)==options.ROIsize && size(corrfn,2)==options.ROIsize % second two arguments added MM because of error when corrfn size is strange
                %do plain symmetric gaussian fit to data
                [coeffGtime] = gaussfit(corrfn,'time',options.pixelSize,'n',options.fitRadius);
                [coeffGrotated] = gaussfit(corrfn,'rotated',options.pixelSize,'n',options.fitRadius);
                velocityMap{k}.coeffGtime{i,j} = coeffGtime;
                velocityMap{k}.coeffGrotated{i,j} = coeffGrotated;
                % Discards some lags of the correlation functions if they have
                % fitted omegas larger than a threshold
                %cutOff =max([1 min([find(coeffGtime(:,2)>options.omegaThreshold,1,'first') find(coeffGtime(:,3)>options.omegaThreshold,1,'first')]) ]);
%                 cutOff =max(1,min([find(coeffGtime(:,2)>options.omegaThreshold,1,'first') find(coeffGtime(:,3)>options.omegaThreshold,1,'first')]));
                
                % set beam waists above cutoff to zero so you can
                % see where it got cut off
%                 coeffGtime(cutOff:end,2:3) = 0;
                % crop time corr as well
                %corrfn(:,:,cutOff:end) = 0;
                [velocityMap{k}.Px(u,2:-1:1), velocityMap{k}.Sx(u)] = polyfit(options.timeFrame*(1:size(corrfn,3)), squeeze(coeffGtime(1:size(corrfn,3),5))', 1);
                [velocityMap{k}.Py(u,2:-1:1), velocityMap{k}.Sy(u)] = polyfit(options.timeFrame*(1:size(corrfn,3)), squeeze(coeffGtime(1:size(corrfn,3),6))', 1);
            else
                velocityMap{k}.Px(u,2:-1:1) = NaN; 
%                 velocityMap{k}.Sx(u,:) = NaN;
                
                velocityMap{k}.Py(u,2:-1:1) = NaN; 
%                 velocityMap{k}.Sy(u,:) = NaN;
            end

    end % end for indx

    % remove bad vectors for each vector map
    posx=squeeze(position_x); % What is this for? Unnecesary?
    posx=posx(1:size(posx,1)*size(posx,2));
    posy=squeeze(position_y);
    posy=posy(1:size(posy,1)*size(posy,2));
    vx=nan(size(position_x));
    vy=nan(size(position_y));
    for u=1:length(iroi)
        i=iroi(u);
        j=jroi(u);
    vy(i,j)=squeeze(velocityMap{k}.Py(u,2));
    vx(i,j)=squeeze(velocityMap{k}.Px(u,2));
    end
    velocityMap{k}.vx=vx;
    velocityMap{k}.vy=vy;
    vx=vx(1:size(vx,1)*size(vx,2));
    vy=vy(1:size(vy,1)*size(vy,2));

    goodVectorsx = vectorOutlier(posx', posy', vx, options.ROIshift, options.threshVector);
    goodVectorsy = vectorOutlier(posx', posy', vy, options.ROIshift, options.threshVector);
    
    velocityMap{k}.goodVectors = (goodVectorsx & goodVectorsy)' & (~isnan(vx) & ~isnan(vy)); 
%     clear goodVectorsx goodVectorsy
    % Throw out vectors whose magnitude is greater than 3 std devs
    % away from the mean
%     for m=1:2 % run twice!
%         magnitudesPerSec = sqrt(vx(velocityMap{k}.goodVectors).^2+vy(velocityMap{k}.goodVectors).^2)';
%         badVectors = find(abs(sqrt(vx.^2+vy.^2)-mean(magnitudesPerSec))>(std(magnitudesPerSec)*3));
%         velocityMap{k}.goodVectors(badVectors) = 0; %#ok<FNDSB>
%     end
    
    velocityMap{k}.goodVectors = reshape(velocityMap{k}.goodVectors, along_x, along_y);
    velocityMap{k}.Vx = reshape(velocityMap{k}.Px(:,2), [along_x, along_y, 2]);
    velocityMap{k}.Vy = reshape(velocityMap{k}.Py(:,2), [along_x, along_y, 2]);
    velocityMap{k}.VelMap = sqrt(velocityMap{k}.vx.^2 + velocityMap{k}.vy.^2);
    
    tEnd = toc(tStart);
    fprintf(' ____ TOI finished in %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60)); % added MM
    
end


% output data to file
save([opt.path 'VelocityMap' opt.outputName '.mat'],'velocityMap',  '-v7.3');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_x','-append', '-v7.3');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_y','-append', '-v7.3');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_t','-append', '-v7.3');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'opt','-append', '-v7.3');

end 

