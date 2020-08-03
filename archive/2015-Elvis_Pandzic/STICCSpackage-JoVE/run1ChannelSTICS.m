%% LOAD DATA INTO WORKSPACE
% load the file you would like to analyze
%%
 [filename, pathname] = uigetfile('*','Please select your data set');
opt.pathin=pathname;
ch1name=filename(1:end-4);
% if data is in MultiTiff files
if strcmp(filename(end-2:end),'tif')
image_data1=rd_img16([opt.pathin ch1name '.tif']);
croped1=image_data1;
% if data is in microscope formats such as Zeiss .lsm files, 'imreadBF' to
% load each channel data. Please open the 'imreadBF.m' and see how the
% input parameters are defined in order to load your data. 
elseif strcmp(filename(end-2:end),'lsm')
image_data1=imreadBF([opt.pathin ch1name '.lsm'],1,1:100,1);
croped1=image_data1;
% if data is stored in a .mat variable
elseif strcmp(filename(end-2:end),'mat')
image_data1=load([opt.pathin ch1name '.mat']);
croped1=image_data1;
else 
   error('Your input data set should be either ''.tif'' ,''.lsm'' or ''.mat'' '); 
end
%% OPTIONAL: IF NECESSARY, CROP SMALLER RECTANGLE IN YOUR IMAGE 
%%% CONTAINING ONLY THE AREA OF INTEREST
[croped1,rect] = serimcropold(image_data1);
%% OPTIONAL: DISPLAY YOUR IMAGE SERIES USING THE STACK VIEWER
sv(croped1,'c',5)
%% run STICS
% define a new directory in which STICS outputs will be saved
mkdir([opt.pathin 'STICSOutputs'])
opt.path =[opt.pathin 'STICSOutputs\'];
[velocityMap,position_x, position_y,position_t,opt]=stics_vectormapping(croped1,opt);
%% LOADING STICS RESULTS (if they are not already in workspace)
opt.pathin='C:\Users\George\Desktop\STICCSpackage\';
opt.path =[opt.pathin 'STICSOutputs\'];
load([opt.path 'VelocityMapF-actin Flow']);
%% DISPLAY VECTOR MAPS
% This defines at what velocity um/min do you want to threshold your vector
% maps...set = Inf, if you do not want to threshold. Thresholding will set
% all the ROI with v>threshold to a 'bad' vector and will not plot it...
% The suggested value would be the theoretically maximal possible flow
% speed that one can extract for a given ROI size, pixels size, frame time
opt.thresholdV=((sqrt((opt.ROIsize^2)/2))*opt.pixelSize/(opt.timeFrame*opt.tauLimit))*60; %um/min
opt.timerDisplay='LL' % option for placing the timer in the movie/images
%'LR'=lower right corner, 'UR'=upper right, 'LL'=lower left corner and
%'UL'=upperleft...if set to anything else than no timer is displayed in
%output maps and movie
% Here specify how often do you want to save the images of maps...in case
% you do not want to save them all 
opt.OutputEvery=50;
velocityMap=plotSingleVectorMapOnImage(velocityMap,position_t,position_x,position_y,opt);
%% PLOT VELOCITY MAGNITUDE HISTOGRAM 
opt.thresholdV=Inf;
dummy=0;
for i=1:length(position_t)
  vy=squeeze(velocityMap{i}.vy);
  vx=squeeze(velocityMap{i}.vx);
  velocityMap{1,i}.magnitudesPerMin=sqrt(vx.^2+vy.^2)*60;
velmag=velocityMap{1,i}.magnitudesPerMin;
numvel=size(velmag,1)*size(velmag,2);
velMag(1+dummy:dummy+numvel)=velmag;
dummy=dummy+numvel;
end
hist(log10(velMag),length(velMag)/50)
%set(gcf,'Color ',[1 1 1])
%set(gca,'FontSize',25)
xlim([-2 2])
xlabel('log_{10} |v| (\mum / min)','FontSize',25)
ylabel('Occurence','FontSize',25)
export_fig([opt.path 'VelMagHistogram' opt.outputName '.pdf'])
close
% hist((velMag),length(velMag)/50)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',25)
% %xlim([-2 2])
% xlabel('|v| (\mum / min)','FontSize',25)
% ylabel('Occurence','FontSize',25)
% export_fig([opt.path 'VelMagHistogramNonLoged' opt.outputName '.pdf'])
% close

%%%%%%%%%%%%%% BATCH PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%
% if you have many files to process, you can run the batch process. First you will need to 
% select the polygons containing the regions of interest for each time series, which are 
% saved, and then the STICS is run sequentially on every file inside the given folder.
%%
% First define all the files in the folder that you want the process. If all files are '.tif' 
%you can automatically define the list as follows:
list=dir([opt.pathin '*.tif']);
% if the files are '.lsm' than use list=dir([opt.pathin '*.lsm']);
% now use a 'for' loop to load all the files, one-by-one, and select the
% polygon of interest for each series:
for i=1:length(list)
    namefile=list(i).name
image_data=rd_img16([opt.pathin namefile]); 
[postouse,maskCell]=DetermineCellMaskFromAverageImage(opt.pathin,namefile(1:length(namefile)-4),image_data,16,4);
end
%% Proceed with STICS batch processing
FilterString = {'FourierWhole','none'}; 
mkdir([opt.pathin 'STICSBatch'])
opt.path =[opt.pathin 'STICSBatch/'];
for i=1:length(list)
namefile=list(i).name
image_data=rd_img16([opt.pathin namefile]);
opt.filtering = FilterString{i}; 
opt.outputName = namefile(1:length(namefile)-4); 
load([opt.pathin 'MaskWithinCell_Serie' namefile(1:length(namefile)-4) '.mat']);
opt.postouse=postouse;
opt.maskCell=maskCell;
[velocityMap,position_x, position_y,position_t,opt]=stics_batch(image_data,opt);
end

