function [velocityMap,position_x, position_y,position_t,opt]=stics_vectormappingpar(series,opt)
%this functions takes the timeseries and calculates 
%for defined regions of
%size 'ROIsize' and 'TOIsize' stics correlation function 
%and fits symmetric
%Gaussian to extract the flow in x and y directions...
%returns velocity maps, x, y and t positions and opt parameters

%Developed by Elvis Pandzic 2012-2014 under supervision of 
%Prof. Paul W. Wiseman @ McGill University. 
% correspondence addressed to : e.pandzic@unsw.edu.au


% Look at opt structure and either assign default value, or use passed
% value
if nargin==1
    opt.pixelSize = 0.1; % just initialize it 
    %to something so the code below doesn't throw an error
end

%these parameters can be defined in the run1channel before calling vectmap
%not used to see ifparticular field already defined or not
if ~isfield(opt, 'pixelSize'), opt.pixelSize = 0.16; end 
%size of pixel in  microns
if ~isfield(opt, 'timeFrame'), opt.timeFrame = 0.3; end 
%time delay (s) between subsequent frames
if ~isfield(opt, 'tauLimit'), opt.tauLimit = 6; end
% what is highest time lag to compute stics of TO
%maximum number of time lags upto which time correlation calculated
%#Frames/5
%filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
if ~isfield(opt, 'filtering'), opt.filtering = 'FourierWhole'; end
%if ~isfield(opt, 'MoveAverage'), opt.MoveAverage = 21; end
%removing spurious vectors
if ~isfield(opt, 'fitRadius'), opt.fitRadius = 8; end
%how many pixels (radius) are considered around the peak of 
%corr fn when fitting 2D Gaussian. 
% radius of circular mask centered at local maxima in ROI.
if ~isfield(opt, 'omegaThreshold'), opt.omegaThreshold = 50; end
% % how big (pixels) do you allow the radius of corr fn 
%to be...removes wide corr fns...if much larger than diff limited PSF
if ~isfield(opt, 'threshVector'), opt.threshVector = 2; end 
% what is the threshold delta(v) (um/s) used in 
%discarding spurious vectors. significant overlap between data analysed
%in neighbouring regions. discarded if the value exceeds by this amount
%from neighbour
%if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
if ~isfield(opt, 'ROIsize'), opt.ROIsize = 16; end
% ROI size in pixels
%Each arrow on the final velocity map is calculated from a 
%region of pixels (ROIsize)x(ROIsize).
if ~isfield(opt, 'ROIshift'), opt.ROIshift = 4; end
%what is ROI centers shift in pixels
if ~isfield(opt, 'TOIsize'), opt.TOIsize = 5; end 
% what is TOI size in frames;time equilvalent of ROI
%look for change in Image series
if ~isfield(opt, 'TOIshift'), opt.TOIshift = 1; end 
%how much is TOI shifted
if ~isfield(opt, 'axisTitle'), opt.axisTitle = 'Sub30_movie011_HistogramMatched_timecroppedMSD'; end %what will be used on vector map title
if ~isfield(opt, 'outputName'), opt.outputName = 'Sub30_movie011_HistogramMatched_timecroppedMSD'; end %output file name
if ~isfield(opt, 'path'), opt.path ='D:\Pam\MATLAB\STICS\sample imaging data\'; end
if ~isfield(opt, 'exportimages'), opt.exportimages ='y'; end %'y' if you want export a pdf of every vector map, 'n' if you do not
if ~isfield(opt, 'imagesformat'), opt.imagesformat ='png'; end % 'png' or 'pdf' images...can add other option (formats) in plotSingleVectorMap code
if ~isfield(opt, 'movieformat'), opt.movieformat ='mp4'; end % movie format can be avi, jpeg or mp4


% shift of ROI's amount & define the position of ROI-TOI
fracROIshift=opt.ROIsize/opt.ROIshift;
fracTOIshift=opt.TOIsize/opt.TOIshift;

%how many roi/toi in the cropped image
%+1 to add the first set of pixels
along_y = floor((size(series,2)-opt.ROIsize)/opt.ROIshift)+1;
along_x = floor((size(series,1)-opt.ROIsize)/opt.ROIshift)+1;
along_t = floor((size(series,3)-opt.TOIsize)/opt.TOIshift)+1;

position_x  =   zeros(1,along_x);
position_y  =   zeros(1,along_y);
position_t  =   zeros(1,along_t);
% defining all of the ROI centers positons in x and y
parfor i=1:along_x
position_x(i) = 1/2 + (i-1)/fracROIshift*opt.ROIsize + opt.ROIsize/2;
end
%duplicates x_position along_y times and then transposes 
position_x=repmat(position_x,along_y,1);
position_x=position_x';
parfor j=1:along_y
position_y(j) = 1/2 + (j-1)/fracROIshift*opt.ROIsize + opt.ROIsize/2;
end
%duplicates y_position along_x times.
%dimensions of x_poistion, y_position same (along_x x along_y)
position_y=repmat(position_y,along_x,1);

%pre-process to determine if you want to analyse the whole FOV or select a
%polygon to analyze...in case of polygon, it will speed up the whole
%calculations because less ROI-TOI to consider...
%imagesc displays image with scaled colours, mean image in time
imagesc(mean(series,3)) 
prompt = {'want to select the polygon or use whole FOV?'};
    dlg_title = 'y or n';
    num_lines = 1;
    def = {'y'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
if strcmp(answer,'y')
maskCell=[];
%roipoly is an inbuilt function in matlab in Image processing toolbox
%x,y define the image limit, poly1 has polygon vertices coordinates
%image inside poly saved in maskCell in binary format
%dimension of mask cell is same as series(cropped image)
%building a polygon mask in the whole image
while isempty(maskCell)
    [x,y,maskCell,poly1.x,poly1.y]=roipoly;
end
elseif strcmp(answer,'n')
    maskCell=ones(size(series,1),size(series,2));
end
close

%selecting ROI size sub arrays from maskCell
postouse=zeros(along_x,along_y);
parfor w=1:along_x
for s=1:along_y
  %checking where the pixels inside polygon present ROIxROI at a time
  indx=1+(w-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((w-1)/fracROIshift)));
  indy=1+(s-1)/fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((s-1)/fracROIshift)));         
  if(maskCell(indx,indy)==1)
  postouse(w,s)=1;
  end
end
end
opt.postouse=postouse;
opt.maskCell=maskCell;
%immobile filter data if necessary
if strcmp(opt.filtering,'FourierWhole')
    series= immfilter(series,'F',1); 
elseif strcmp(opt.filtering,'MovingAverage')
    series= immfilter(series,opt.MoveAverage); 
elseif strcmp(opt.filtering,'butterIIR')
    %set the cutoff frequency to be the 2*sampling frequency / size of the series 
    % this way we will certainly remove the immobile component
    series= butterIIR(series, 2/(opt.timeFrame*size(series,3)),'none');
elseif strcmp(opt.filtering,'none')
else
    error('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
end

% %indexes of ROI within the mask of interest (polygon)
[iroi,jroi]=find(postouse);

for k=1:along_t % loop along time
    % position of TOI centers in time
    %centre time frame
  
    position_t(k) = 1/2 + (k-1)/fracTOIshift*opt.TOIsize + opt.TOIsize/2;
    % define whole FOV per TOI
    
    TOIFOV=series(:,:,floor(1+(k-1)/fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/fracTOIshift))));
    % define the average FOV TOI image 
    %that will be used in the display of the vector maps
    velocityMap{k}.data_raw=mean(TOIFOV,3);
    
    tStart = tic; % 
    fprintf('analyzing TOI %i of %i ____ date and time is %s',k,along_t,datestr(now)); 
    velocityMap_Px1=zeros(length(iroi),2);
    velocityMap_Py1=zeros(length(iroi),2);
    velocityMap_Px2=zeros(length(iroi),2);
    velocityMap_Py2=zeros(length(iroi),2);
 
     parfor u=1:length(iroi) %accessing the data of only 
     %the ROI that are defined within the selected polygon
         i=iroi(u);
         j=jroi(u);
          
            %define region (x,y,t) on which stics will be applied
            %defined by postouse
            regionanalyse=TOIFOV((1+(i-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((i-1)/fracROIshift))),(1+(j-1)/fracROIshift*opt.ROIsize):(opt.ROIsize*(1+((j-1)/fracROIshift))),:);
           
            %apply regular cross-corr stics
            [corrfn]=stics(regionanalyse,opt.tauLimit);
            %[i j k]
         
            %fit vxtau and vytau vs tau linearly to extract vx and vy
            if size(corrfn,3)>1 && size(corrfn,1)==opt.ROIsize && size(corrfn,2)==opt.ROIsize % second two arguments added MM because of error when corrfn size is strange
                %do plain symmetric gaussian fit to data
                %coeffGtime-coefficients from Gaussian fitting 
                [coeffGtime] = gaussfit(corrfn,'time',opt.pixelSize,'n',opt.fitRadius);
                %velocityMap_coeffGtime{k}.coeffGtime=coeffGtime;
                % Discards some lags of the correlation functions if they have
                % fitted omegas larger than a threshold
                %cutOff =max([1 min([find(coeffGtime(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime(:,3)>opt.omegaThreshold,1,'first')]) ]);
                cutOff =max(1,min([find(coeffGtime(:,2)>opt.omegaThreshold,1,'first') find(coeffGtime(:,3)>opt.omegaThreshold,1,'first')]));
                
                % set beam waists above cutoff to zero so you can
                % see where it got cut off
                coeffGtime(cutOff:end,2:3) = 0;
                % crop time corr as well
                %corrfn(:,:,cutOff:end) = 0;
              
                fit1 = polyfit(opt.timeFrame*(1:size(corrfn,3)),squeeze(coeffGtime(1:size(corrfn,3),5))',1);
                fit2 = polyfit(opt.timeFrame*(1:size(corrfn,3)),squeeze(coeffGtime(1:size(corrfn,3),6))',1);
                
                velocityMap_Px1(u) = fit1(1);
                velocityMap_Px2(u)=fit1(2);
                velocityMap_Py1(u) = fit2(1);
                velocityMap_Py2(u) = fit2(2);
               
                
            else
                
                velocityMap_Px1(u) = NaN;
                velocityMap_Px2(u) = NaN;
                velocityMap_Py1(u) = NaN;
                velocityMap_Py2(u) = NaN;
                
            end

      end % end for indx
    
    
     
      velocityMap{k}.Px(:,2) = velocityMap_Px1(:,1);
      velocityMap{k}.Px(:,1) = velocityMap_Px2(:,1);
      velocityMap{k}.Py(:,2) = velocityMap_Py1(:,1);
      velocityMap{k}.Py(:,1) = velocityMap_Py2(:,1);
   
    %velocityMap{k}.coeffGtime=velocityMap_coeffGtime;
    
    
    
    % remove bad vectors for each vector map
    posx=squeeze(position_x);
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
    vy=vy(1:size(vy,1)*size(vy,2));
    vx=vx(1:size(vx,1)*size(vx,2));
    goodVectorsx = vectorOutlier(posx',posy',vx,opt.ROIshift,opt.threshVector);
    goodVectorsy = vectorOutlier(posx',posy',vy,opt.ROIshift,opt.threshVector);
    
    velocityMap{k}.goodVectors = (goodVectorsx & goodVectorsy)' & (~isnan(vx) & ~isnan(vy)); clear goodVectorsx goodVectorsy
    % Throw out vectors whose magnitude is greater than 3 std devs
    % away from the mean
    for m=1:2 % run twice!
        magnitudesPerSec = sqrt(vx(velocityMap{k}.goodVectors).^2+vy(velocityMap{k}.goodVectors).^2);
        badVectors = find(abs(sqrt(vx.^2+vy.^2)-mean(magnitudesPerSec))>(std(magnitudesPerSec)*3));
        velocityMap{k}.goodVectors(badVectors) = 0; %#ok<FNDSB>
    end
    
    
    tEnd = toc(tStart);
    fprintf(' ____ TOI finished in %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60)); % added MM
    
end


% output data to file
save([opt.path 'VelocityMap' opt.outputName '.mat'],'velocityMap');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_x','-append');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_y','-append');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_t','-append');
save([opt.path 'VelocityMap' opt.outputName '.mat'],'opt','-append');

end 

