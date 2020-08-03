function [velocityMap,position_x, position_y,position_t,opt]=sticcs_vectormapping(series1,series2,opt)

%Description: this functions takes the timeseries and calculates for defined regions of
%size 'ROIsize' and 'TOIsize' sticcs (auto- and cross-correlation between two channels) 
%correlation funstion and fits symmetric Gaussian to extract the flow in x and y directions...

%Developed by Elvis Pandzic 2012-2014 under supervision of 
%Prof. Paul W. Wiseman @ McGill University. 
% correspondence addressed to : e.pandzic@unsw.edu.au


tic
% Look at opt structure and either assign default value, or use passed
% value
if nargin==2 
    opt.pixelSize = 0.1; % just initialize it to something so the code below doesn't throw an error
end


if ~isfield(opt, 'pixelSize'), opt.pixelSize = 0.1031 ; end % in micrometers
if ~isfield(opt, 'timeFrame'), opt.timeFrame = 10; end %time delay (s) between subsequent frames
if ~isfield(opt, 'timeDelay'), opt.timeDelay = 0; end %extra time delay (s) between two channels
if ~isfield(opt, 'tauLimit'), opt.tauLimit = 20; end % what is highest time lag to comute stics of TOI
%filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
if ~isfield(opt, 'filtering'), opt.filtering ='FourierWhole'; end 
if ~isfield(opt, 'MoveAverage'), opt.MoveAverage = 21; end
if ~isfield(opt, 'fitRadius'), opt.fitRadius = 5; end  % how many pixels (radius) are considered around the peak of corr fn when fitting 2D Gaussian 
if ~isfield(opt, 'omegaThreshold'), opt.omegaThreshold = 15; end % how big (pixels) do you allow the radius of corr fn to be...removes wide corr fns...
if ~isfield(opt, 'threshVector'), opt.threshVector = 8; end % what is the threshold delta(v) (um/s) used in discarding spurious vectors
%if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
if ~isfield(opt, 'ROIsize'), opt.ROIsize = 16; end % ROI size in pixels
if ~isfield(opt, 'ROIshift'), opt.ROIshift = 4; end %what is ROI centers shift in pixels
if ~isfield(opt, 'TOIsize'), opt.TOIsize = 60; end % what is TOI size in frames
if ~isfield(opt, 'TOIshift'), opt.TOIshift = 1; end %how much is TOI shifted
if ~isfield(opt, 'axisTitle'), opt.axisTitle = 'TwoColorVelocityMap F-actin and PM'; end %what will be used on vector map title
if ~isfield(opt, 'outputName'), opt.outputName = 'TwoColorVelocityMap F-actin and PM'; end %output file name
if ~isfield(opt, 'path'), opt.path ='C:\Users\YOGA-PW\Dropbox (Wiseman Research)\Data_Elvis\STICCSOutput'; end
if ~isfield(opt, 'exportimages'), opt.exportimages ='y'; end %'y' if you want export a pdf of every vector map, 'n' if you do not
if ~isfield(opt, 'imagesformat'), opt.imagesformat ='png'; end % 'png' or 'pdf' images...can add other option (formats) in plotSingleVectorMap code
if ~isfield(opt, 'movieformat'), opt.movieformat ='mp4'; end % movie format can be avi, jpeg or mp4




along_y = floor((size(series1,2)-opt.ROIsize)/opt.ROIshift)+1;
along_x = floor((size(series1,1)-opt.ROIsize)/opt.ROIshift)+1;
along_t = floor((size(series1,3)-opt.TOIsize)/opt.TOIshift)+1;

% shift of ROI's amount
fracROIshift=opt.ROIsize/opt.ROIshift;
fracTOIshift=opt.TOIsize/opt.TOIshift;
position_x  =   zeros(1,along_x);
position_y  =   zeros(1,along_y);
position_t  =   zeros(1,along_t);

% defining all of the ROI centers positons in x and y
for i=1:along_x
position_x(i) = 1/2 + (i-1)/fracROIshift*opt.ROIsize +  opt.ROIsize/2;
end
position_x=repmat(position_x,along_y,1);
position_x=position_x';
for j=1:along_y
position_y(j) = 1/2 + (j-1)/fracROIshift*opt.ROIsize + opt.ROIsize/2;
end
position_y=repmat(position_y,along_x,1);



%pre-process to determine if you want to analyse the whole FOV or select a
%polygon to analyze...in case of polygon, it will speed up the whole
%calcualations because less ROI-TOI to consider...
imcomp(mean(series1,3),mean(series2,3),'y'); 
prompt = {'want to select the polygon or use whole FOV?'};
    dlg_title = 'y or n';
    num_lines = 1;
    def = {'y'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
if strcmp(answer,'y')
maskCell=[];
while isempty(maskCell)
    [x,y,maskCell,poly1.x,poly1.y]=roipoly;
end
elseif strcmp(answer,'n')
    maskCell=ones(size(series1,1),size(series1,2));
end
close
postouse=zeros(along_x,along_y);
for w=1:along_x
for s=1:along_y
  indx=1+(w-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((w-1)/fracROIshift)));
  indy=1+(s-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((s-1)/fracROIshift)));  
  indx = floor(indx);
  indy = floor(indy);
if(maskCell(indx,indy)==1)
postouse(w,s)=1;
end
end
end
opt.postouse=postouse;
opt.maskCell=maskCell;

%immobile filter data
if strcmp(opt.filtering,'FourierWhole')
   series1= immfilter(series1,'F',1);
   series2= immfilter(series2,'F',1);
elseif strcmp(opt.filtering,'MovingAverage')
   series1= immfilter(series1,opt.MoveAverage);
   series2= immfilter(series2,opt.MoveAverage);
elseif strcmp(opt.filtering,'butterIIR')
    %set the cutoff frequency to be the 2*sampling frequency / size of the series 
    % this way we will certainly removed the immobile component
   series1= butterIIR(series1, 2/(opt.timeFrame*size(series1,3)),'none');
   series2= butterIIR(series2, 2/(opt.timeFrame*size(series2,3)),'none');
elseif strcmp(opt.filtering,'none')
else
       error('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
         
end

% %indexes of ROI within the mask of interest (polygon)
[iroi jroi]=find(postouse);

for k=1:along_t % loop along time
    % position of TOI centers in time
    position_t(k) = 1/2 + (k-1)/fracTOIshift*opt.TOIsize + opt.TOIsize/2;
      % define whole FOV TOI
    TOIFOV1=series1(:,:,floor(1+(k-1)/fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/fracTOIshift))));
    TOIFOV2=series2(:,:,floor(1+(k-1)/fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/fracTOIshift)))); 
    % define the average FOV TOI image that will be used in the display of the vector maps
    velocityMap{1,k}.data_raw=mean(TOIFOV1,3);
    velocityMap{2,k}.data_raw=mean(TOIFOV2,3);
  
    
    tStart = tic; % added MM
    fprintf('analyzing TOI %i of %i ____ date and time is %s',k,along_t,datestr(now)); 
   
    
  for u=1:length(iroi) %accessing the data of only 
     %the ROI that are defined within the selected polygon
         i=iroi(u);
         j=jroi(u);  
         %define regions (x,y,t) on which sticcs will be applied  
        regionanalyse1=TOIFOV1((1+(i-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((i-1)/fracROIshift))),(1+(j-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((j-1)/fracROIshift))),:);
        regionanalyse2=TOIFOV2((1+(i-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((i-1)/fracROIshift))),(1+(j-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((j-1)/fracROIshift))),:);
                  
        %apply regular sticcs 
        [corrfn12 corrfn21 corrfn1 corrfn2]=stics(regionanalyse1,regionanalyse2,opt.tauLimit);
         
        if size(corrfn1,3)>1
        %do plain symmetric gaussian fit to data
        [coeffGtime1] = gaussfit(corrfn1,'time',opt.pixelSize,'n',opt.fitRadius);
        % Discards some lags of the correlation functions if they have
        % fitted omegas larger than a threshold
        cutOff1 =max(1,min([find(coeffGtime1(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime1(:,3)>opt.omegaThreshold,1,'first') ]));
        % set beam waists above cutoff to zero so you can
        % see where it got cut off
        coeffGtime1(cutOff1:end,2:3) = 0;
        velocityMap{1,k}.Px(u,2:-1:1) = polyfit((opt.timeFrame+opt.timeDelay)*(0:size(corrfn1,3)-1),squeeze(coeffGtime1(1:size(corrfn1,3),5))',1);
        velocityMap{1,k}.Py(u,2:-1:1) = polyfit((opt.timeFrame+opt.timeDelay)*(0:size(corrfn1,3)-1),squeeze(coeffGtime1(1:size(corrfn1,3),6))',1);
        else
         velocityMap{1,k}.Px(u,1:2) =NaN;
         velocityMap{1,k}.Py(u,1:2) =NaN;
        end
        % do it for channel 2
        if size(corrfn2,3)>1
        %do plain symmetric gaussian fit to data
        [coeffGtime2] = gaussfit(corrfn2,'time',opt.pixelSize,'n',opt.fitRadius);
        % Discards some lags of the correlation functions if they have
        % fitted omegas larger than a threshold
        cutOff2 =max(1,min([find(coeffGtime2(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime2(:,3)>opt.omegaThreshold,1,'first') ]));
        % set beam waists above cutoff to zero so you can
        % see where it got cut off
        coeffGtime2(cutOff2:end,2:3) = 0;
        velocityMap{2,k}.Px(u,2:-1:1) = polyfit((opt.timeFrame+opt.timeDelay)*(0:size(corrfn2,3)-1),squeeze(coeffGtime2(1:size(corrfn2,3),5))',1);
        velocityMap{2,k}.Py(u,2:-1:1) = polyfit((opt.timeFrame+opt.timeDelay)*(0:size(corrfn2,3)-1),squeeze(coeffGtime2(1:size(corrfn2,3),6))',1);
        else
         velocityMap{2,k}.Px(u,1:2) =NaN;
         velocityMap{2,k}.Py(u,1:2) =NaN;
        end
        % do it for cross-corr 12
        if size(corrfn12,3)>1
        %do plain symmetric gaussian fit to data
        [coeffGtime12] = gaussfit(corrfn12,'time',opt.pixelSize,'n',opt.fitRadius);
         % Discards some lags of the correlation functions if they have
        % fitted omegas larger than a threshold
        cutOff12 =max(1,min([find(coeffGtime12(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime12(:,3)>opt.omegaThreshold,1,'first') ]));
        % set beam waists above cutoff to zero so you can
        % see where it got cut off
        coeffGtime12(cutOff12:end,2:3) = 0;
        velocityMap{3,k}.Px(u,2:-1:1) = polyfit(opt.timeDelay+((opt.timeFrame+opt.timeDelay)*(0:size(corrfn12,3)-1)),squeeze(coeffGtime12(1:size(corrfn12,3),5))',1);
        velocityMap{3,k}.Py(u,2:-1:1) = polyfit(opt.timeDelay+((opt.timeFrame+opt.timeDelay)*(0:size(corrfn12,3)-1)),squeeze(coeffGtime12(1:size(corrfn12,3),6))',1);
        else
         velocityMap{3,k}.Px(u,1:2) =NaN;
         velocityMap{3,k}.Py(u,1:2) =NaN;
        end
        % do it for cross-corr 21
        if size(corrfn21,3)>1
        %do plain symmetric gaussian fit to data
        [coeffGtime21] = gaussfit(corrfn21,'time',opt.pixelSize,'n',opt.fitRadius);
        % Discards some lags of the correlation functions if they have
        % fitted omegas larger than a threshold
        cutOff21 =max(1,min([find(coeffGtime21(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime21(:,3)>opt.omegaThreshold,1,'first') ]));
        % set beam waists above cutoff to zero so you can
        % see where it got cut off
        coeffGtime21(cutOff21:end,2:3) = 0;
        velocityMap{4,k}.Px(u,2:-1:1) = polyfit(opt.timeDelay+((opt.timeFrame+opt.timeDelay)*(0:size(corrfn21,3)-1)),squeeze(coeffGtime21(1:size(corrfn21,3),5))',1);
        velocityMap{4,k}.Py(u,2:-1:1) = polyfit(opt.timeDelay+((opt.timeFrame+opt.timeDelay)*(0:size(corrfn21,3)-1)),squeeze(coeffGtime21(1:size(corrfn21,3),6))',1);
        else
         velocityMap{4,k}.Px(u,1:2) =NaN;
         velocityMap{4,k}.Py(u,1:2) =NaN;
        end
        
         
  end % end 'u' loop
  
  
% % % remove bad vectors for each channel and cross-corrs
  for m=1:4 
    posx=squeeze(position_x);
    posx=posx(1:size(posx,1)*size(posx,2));
    posy=squeeze(position_y);
    posy=posy(1:size(posy,1)*size(posy,2));
    vx=nan(size(position_x));
    vy=nan(size(position_y));
    for u=1:length(iroi)
        i=iroi(u);
        j=jroi(u);
    vy(i,j)=squeeze(velocityMap{m,k}.Py(u,2));
    vx(i,j)=squeeze(velocityMap{m,k}.Px(u,2));
    end
    velocityMap{m,k}.vx=vx;
    velocityMap{m,k}.vy=vy;
    vy=vy(1:size(vy,1)*size(vy,2));
    vx=vx(1:size(vx,1)*size(vx,2));
    goodVectorsx = vectorOutlier(posx',posy',vx,opt.ROIshift,opt.threshVector);
    goodVectorsy = vectorOutlier(posx',posy',vy,opt.ROIshift,opt.threshVector);
    velocityMap{m,k}.goodVectors = (goodVectorsx & goodVectorsy)' & (~isnan(vx) & ~isnan(vy)); clear goodVectorsx goodVectorsy
    % Throw out vectors whose magnitude is greater than 3 std devs
    % away from the mean of the neigbouring vectors
    for j=1:2 % run twice!
    magnitudesPerSec = sqrt(vx(velocityMap{m,k}.goodVectors).^2+vy(velocityMap{m,k}.goodVectors).^2);
    badVectors = find(abs(sqrt(vx.^2+vy.^2)-mean(magnitudesPerSec))>(std(magnitudesPerSec)*3));
    velocityMap{m,k}.goodVectors(badVectors) = 0; %#ok<FNDSB>
    end
  end %end 'm' loop
  
  
  tEnd = toc(tStart); 
  fprintf(' ____ TOI finished in %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60)); 
    
  
end %end 'k' loop

% output data to file
save([opt.path 'TwoColorVelocityMap' opt.outputName '.mat'],'velocityMap');
save([opt.path 'TwoColorVelocityMap' opt.outputName '.mat'],'position_x','-append');
save([opt.path 'TwoColorVelocityMap' opt.outputName '.mat'],'position_y','-append');
save([opt.path 'TwoColorVelocityMap' opt.outputName '.mat'],'position_t','-append');
save([opt.path 'TwoColorVelocityMap' opt.outputName '.mat'],'opt','-append');


toc

