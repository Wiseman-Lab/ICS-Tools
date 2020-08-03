%% LOAD DATA INTO WORKSPACE
% define files to analyze...hold Shift key while selecting 2 files
[filename, pathname,filterindex] = uigetfile('*','Please select your 2 data sets','MultiSelect', 'on');
opt.pathin=pathname;
ch1name=filename{1}(1:end-4);
ch2name=filename{2}(1:end-4);
% if data is in MultiTiff files
if strcmp(filename{1}(end-2:end),'tif')
image_data1=rd_img16([opt.pathin ch1name '.tif']);
image_data2=rd_img16([opt.pathin ch2name '.tif']);
croped1=image_data1;
croped2=image_data2;
% if data is in microscope formats such as Zeiss .lsm files, 'imreadBF' to
% load each channel data. Please open the 'imreadBF.m' and see how the
% input parameters are defined in order to load your data. 
elseif strcmp(filename{1}(end-2:end),'lsm')
image_data1=imreadBF([opt.pathin ch1name '.lsm'],1,1:100,1);
image_data2=imreadBF([opt.pathin ch2name '.lsm'],1,1:100,1);
croped1=image_data1;
croped2=image_data2;
% if data is stored in a .mat variable
elseif strcmp(filename{1}(end-2:end),'mat')
image_data1=load([opt.pathin ch1name '.mat']);
image_data2=load([opt.pathin ch2name '.mat']);
croped1=image_data1;
croped2=image_data2;
else 
   error('Your input data sets should be either ''.tif'' ,''.lsm'' or ''.mat'' '); 
end
%% OPTIONAL: IF NECESSARY, CROP SMALLER RECTANGLE IN YOUR IMAGE 
%%% CONTAINING ONLY THE AREA OF INTEREST
[croped1,rect] = serimcropold(image_data1);
clear croped2
% crop same rectangle within channel 2
for i=1:size(image_data2,3)
    croped2(:,:,i) = imcrop(image_data2(:,:,i),rect);
end
clear image_data1 image_data2 

%% OPTIONAL: DISPLAY YOUR IMAGE SERIES USING THE STACK VIEWER
sv(croped1,croped2,'c',5)
%% run STICCS
% define a new directory in which STICCS outputs will be saved
mkdir([opt.pathin 'STICCSOutputs'])
opt.path =[opt.pathin 'STICCSOutputs/'];
[velocityMap,position_x, position_y,position_t,opt]=sticcs_vectormapping(croped1,croped2,opt);


%% DISPLAY VECTOR MAPS
%%% This part is run to produce vector maps for 2 color data.
%%% For each auto- and cross-correlation channels, there is a 
%%% separate set of commands below. Make sure to load the output 
%%% of the STICCS before you run lines below 
%% LOADING STICCS RESULTS (if they are not already in workspace)
opt.pathin='C:\Users\George\Desktop\STICCSpackage\';
opt.path =[opt.pathin 'STICCSOutputs\'];
load([opt.path 'TwoColorVelocityMap F-actin and PM']);
%% CHANNEL 1 VECTOR MAP
clear velocityMap2 opt2
opt2=opt;
opt2.axisTitle='Channel 1';
opt2.outputName = ['Channel1_' opt.outputName]; 
opt2.exportimages='y';
% This defines at what velocity um/min do you want to threshold your vector
% maps...set = Inf, if you do not want to threshold. Thresholding will set
% all the ROI with v>threshold to a 'bad' vector and will not plot it...
% The suggested value would be the theoretically maximal possible flow
% speed that one can extract for a given ROI size, pixels size, frame time
opt2.thresholdV=((sqrt((opt.ROIsize^2)/2))*opt.pixelSize/(opt.timeFrame*opt.tauLimit))*60; %um/min
for k=1:length(position_t)
velocityMap2{k}=velocityMap{1,k};
velocityMap2{k}.data_raw=velocityMap{1,k}.data_raw;
end
opt2.timerDisplay='LL' % option for placing the timer in the movie/images
%'LR'=lower right corner, 'UR'=upper right, 'LL'=lower left corner and
%'UL'=upperleft...if set to anything else than no timer is displayed in
%output maps and movie
velocityMap2=plotSingleVectorMapOnImage(velocityMap2,position_t,position_x,position_y,opt2);
%% PLOT VELOCITY MAGNITUDE HISTOGRAM 
dummy=0;
for i=1:length(position_t)
velmag=velocityMap2{1,i}.magnitudesPerMin;
numvel=size(velmag,1)*size(velmag,2);
velMag(1+dummy:dummy+numvel)=velmag;
dummy=dummy+numvel;
end
hist(log10(velMag),length(velMag)/50)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',25)
xlim([-2 0])
xlabel('log_{10} |v| (\mum / min)','FontSize',25)
ylabel('Occurence','FontSize',25)
export_fig([opt2.path 'VelMagHistogram' opt2.outputName '.pdf'])
close

%% CHANNEL 2 VECTOR MAP
clear velocityMap2 opt2
opt2=opt;
opt2.axisTitle='Channel 2';
opt2.outputName = ['Channel2_' opt.outputName];
opt2.exportimages='y';
% This defines at what velocity um/min do you want to threshold your vector
% maps...set = Inf, if you do not want to threshold. thresholding will set
% all the ROI with v>threshold to a 'bad' vector and will not plot it...
% The suggested value would be the theoretically maximal possible flow
% speed that one can extract for a given ROI size, pixels size, frame time
opt2.thresholdV=((sqrt((opt.ROIsize^2)/2))*opt.pixelSize/(opt.timeFrame*opt.tauLimit))*60; %um/min

for k=1:length(position_t)
velocityMap2{k}=velocityMap{2,k};
velocityMap2{k}.data_raw=velocityMap{2,k}.data_raw;
end
opt2.timerDisplay='LL' % option for placing the timer in the movie/images
%'LR'=lower right corner, 'UR'=upper right, 'LL'=lower left corner and
%'UL'=upperleft...if set to anything else than no timer is displayed in
%output maps and movie
velocityMap2=plotSingleVectorMapOnImage(velocityMap2,position_t,position_x,position_y,opt2);

%% PLOT VELOCITY MAGNITUDE HISTOGRAM 
dummy=0;
for i=1:length(position_t)
velmag=velocityMap2{1,i}.magnitudesPerMin;
numvel=size(velmag,1)*size(velmag,2);
velMag(1+dummy:dummy+numvel)=velmag;
dummy=dummy+numvel;
end
hist(log10(velMag),length(velMag)/50)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',25)
xlim([-2 0])
xlabel('log_{10} |v| (\mum / min)','FontSize',25)
ylabel('Occurence','FontSize',25)
export_fig([opt2.path 'VelMagHistogram' opt2.outputName '.pdf'])
close

%% CROSS-CORRELATION 12 (i.e. channel 1 is the reference)
clear velocityMap2 opt2
opt2=opt;
opt2.axisTitle='Cross-Corr-12';
opt2.outputName = ['CC12_' opt.outputName];
opt2.exportimages='y';
% This defines at what velocity um/min do you want to threshold your vector
% maps...set = Inf, if you do not want to threshold. thresholding will set
% all the ROI with v>threshold to a 'bad' vector and will not plot it...
% The suggested value would be the theoretically maximal possible flow
% speed that one can extract for a given ROI size, pixels size, frame time
opt2.thresholdV=((sqrt((opt.ROIsize^2)/2))*opt.pixelSize/(opt.timeFrame*opt.tauLimit))*60; %um/min

opt2.timerDisplay='LL' % option for placing the timer in the movie/images
%'LR'=lower right corner, 'UR'=upper right, 'LL'=lower left corner and
%'UL'=upperleft...if set to anything else than no timer is displayed in
%output maps and movie
for k=1:length(position_t)
velocityMap2{k}=velocityMap{3,k};
velocityMap2{k}.data_raw=(velocityMap{1,k}.data_raw+velocityMap{2,k}.data_raw)/2;
end
velocityMap2=plotSingleVectorMapOnImage(velocityMap2,position_t,position_x,position_y,opt2);


%% PLOT VELOCITY MAGNITUDE HISTOGRAM 
dummy=0;
for i=1:length(position_t)
velmag=velocityMap2{1,i}.magnitudesPerMin;
numvel=size(velmag,1)*size(velmag,2);
velMag(1+dummy:dummy+numvel)=velmag;
dummy=dummy+numvel;
end
hist(log10(velMag),length(velMag)/50)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',25)
xlim([-2 0])
xlabel('log_{10} |v| (\mum / min)','FontSize',25)
ylabel('Occurence','FontSize',25)
export_fig([opt2.path 'VelMagHistogram' opt2.outputName '.pdf'])
close
%% CROSS-CORRELATION 21 (i.e. channel 2 is the reference)
clear velocityMap2 opt2
opt2=opt;
opt2.axisTitle='Cross-Corr-21';
opt2.outputName = ['CC21_' opt.outputName];
opt2.exportimages='y';
% This defines at what velocity um/min do you want to threshold your vector
% maps...set = Inf, if you do not want to threshold. thresholding will set
% all the ROI with v>threshold to a 'bad' vector and will not plot it...
% The suggested value would be the theoretically maximal possible flow
% speed that one can extract for a given ROI size, pixels size, frame time
opt2.thresholdV=((sqrt((opt.ROIsize^2)/2))*opt.pixelSize/(opt.timeFrame*opt.tauLimit))*60; %um/min
opt2.timerDisplay='LL' % option for placing the timer in the movie/images
%'LR'=lower right corner, 'UR'=upper right, 'LL'=lower left corner and
%'UL'=upperleft...if set to anything else than no timer is displayed in
%output maps and movie
for k=1:length(position_t)
velocityMap2{k}=velocityMap{4,k};
velocityMap2{k}.data_raw=(velocityMap{1,k}.data_raw+velocityMap{2,k}.data_raw)/2;
end
velocityMap2=plotSingleVectorMapOnImage(velocityMap2,position_t,position_x,position_y,opt2);

%% PLOT VELOCITY MAGNITUDE HISTOGRAM 
dummy=0;
for i=1:length(position_t)
velmag=velocityMap2{1,i}.magnitudesPerMin;
numvel=size(velmag,1)*size(velmag,2);
velMag(1+dummy:dummy+numvel)=velmag;
dummy=dummy+numvel;
end
hist(log10(velMag),length(velMag)/50)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',25)
xlim([-2 0])
xlabel('log_{10} |v| (\mum / min)','FontSize',25)
ylabel('Occurence','FontSize',25)
export_fig([opt2.path 'VelMagHistogram' opt2.outputName '.pdf'])
close


