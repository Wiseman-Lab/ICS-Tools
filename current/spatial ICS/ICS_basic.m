%% Simple ICS
% To use this code just press 'Run' and then insert the necessary
% information when prompted
% Selection of a file
[file, path] = uigetfile('*.tif','File Selector');
imageSeries = rd_img16([path file]);

%Uncomment this if you would like to just insert the location and file
%name directly
% path = 'D:\ICS feb\serum +- control\';
% imageSeries=rd_img16([path 'cell_40.tif']);

%Uncomment this if you are using the Simulator to generate an image series
% imageSeries = SimulatorImageSeries2();
%%

%Main prompt, most information is gathered here
prompt = {'Start Frame','End Frame','White Noise? (y/n)','Smear Correct? (y/n)','Time Averaging? (y/n)','Pixel Size (um)','Fit Radius (pixels)',...
    'PSF (nm) (G1-344 G1.6-232s R1-373 R1.6-270)','Repetitions','Scanning? (y/n)','Output all values? (or averaged) (y/n)'};
dlg_title = 'Enter Analysis Parameters';
defaultans = {num2str(2),num2str(11),'y','y','y',num2str(0.08125),num2str(12),num2str(344),num2str(3),'n','n'};
answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
sFrame = str2double(answer{1});
fFrame =str2double(answer{2});
wn =answer{3};
smear =answer{4};
timeavg = answer{5};
pixelsize =str2double(answer{6});
fitRadius = str2double(answer{7});
realW = str2double(answer{8})./sqrt(2);
reps = str2double(answer{9});
scan = answer{10};
allV = answer{11};

% Secondary prompt
if strcmp(scan,'y')
    prompt = {'ROI Size (pixels)','Display Graphs? (y/n)'};
    dlg_title = 'Scanning Parameter';
    defaultans = {num2str(64),'n'};
    answer = inputdlg(prompt,dlg_title,[1 60], defaultans);
    ROIsize = str2double(answer{1});
    graphs = answer{2};
else
    prompt = {'Display Graphs? (y/n)'};
    dlg_title = '';
    defaultans = {'y'};
    answer = inputdlg(prompt,dlg_title,[1 60], defaultans);
    graphs = answer{1};
end

% A filter is selected from 4 options, based on the options specific
% filtering of the data will be done
% Note: If you select to time average, then filtering will not be performed
% for individual measurements (there is no point)
filters = {'No Filtering','Outside 3 std','30% around real w','Negative Values'};
[indx,tf] = listdlg('PromptString',{'Select a filter method.','If "no filtering" is selected none of the other ones','will be used.'}...
    ,'ListString',filters,'Name','Filter Selection','ListSize',[250 150],'InitialValue',[2 4]);

%Selection of frames as determined in initial prompt
imageSeries = imageSeries(:,:,sFrame:fFrame);

%White noise correction. Do if you have any background noise
if strcmp(wn,'y')
    [imageSeriesCorrected, noise] = wnCorr(imageSeries);
else
    imageSeriesCorrected = imageSeries;
end
%Unsmearing of TIRF images. Check with current TIRF users if this is
%necessary
if strcmp(smear,'y')
    unsmeared = smearcorr(imageSeriesCorrected, 0.05, 0.0000006); %I don't know where these parameters come from!!
else
    unsmeared = imageSeriesCorrected;
end
%Will average the image series across the third dimension (time)
if strcmp(timeavg,'y')
    unsmeared = mean(unsmeared,3);
end

%Initial crop to select a general area
[croppedImageSeries1, crop] = serimcropold(unsmeared);
%Assignment of structure for final data
if strcmp(scan,'y')
    data = struct('CD_Av',{},'CD_Std',{},'BeamArea_Av',{},'BeamArea_Std',{},'w_Av',{},'w_Std',{},'I_Av',{},'I_Std',{},'R2_Av',{},'R2_Std',{});
else
    data = struct('X_Size',{},'Y_Size',{},'CD_Av',{},'CD_Std',{},'BeamArea_Av',{},'BeamArea_Std',{},'w_Av',{},'w_Std',{},'I_Av',{},'I_Std',{});
end
%%
% croppedImageSeries1 = meanDisp(croppedImageSeries1,0.232,3);

%Will run multiple times allowing you to select multiple ROIs
for j = 1:reps
    %Cropping the image series again
    [croppedImageSeries, crop] = serimcropold(croppedImageSeries1);
    %If scanning ICS is selected this section will run, otherwise proceed to
    %'else'
    if strcmp(scan,'y')
        ICS2DCorrScan = scanICS(croppedImageSeries, crop, fitRadius, pixelsize, ROIsize,realW, wn, graphs, indx,allV);
        %scICSfire(ICS2DCorrScan, croppedImageSeries, crop, 32);
        %This section takes the values from the struct returned in the
        %scanICS code
        particles = zeros(size(ICS2DCorrScan,2),1);
        beamradius = zeros(size(ICS2DCorrScan,2),1);
        beamarea = zeros(size(ICS2DCorrScan,2),1);
        intensities = zeros(size(ICS2DCorrScan,2),1);
        R2s = zeros(size(ICS2DCorrScan,2),1);
        for i = 1:size(particles)
            particles(i) = mean(ICS2DCorrScan(i).CD,'omitnan');
            beamarea(i) = mean(ICS2DCorrScan(i).BeamArea,'omitnan');
            beamradius(i) = mean(ICS2DCorrScan(i).BeamRadius,'omitnan');
            intensity = zeros(size(croppedImageSeries,3),1);
            for j = 1:size(croppedImageSeries,3)
                intensity(j) = mean2(ICS2DCorrScan(i).raw_data(:,:,j));
            end
            intensities(i) = mean(intensity);
            R2s(i) = mean(ICS2DCorrScan(i).R2,'omitnan');
        end
        %If you want all the values then the average for each ROI within
        %the main ROI will be given
        if strcmp(allV,'n')
            data(j).CD_Av = mean(particles,'omitnan');
            data(j).CD_Std = std(particles,'omitnan');
            data(j).BeamArea_Av = mean(beamradius,'omitnan');
            data(j).BeamArea_Std = std(beamarea,'omitnan');
            data(j).w_Av = mean(beamradius,'omitnan');
            data(j).w_Std = std(beamradius,'omitnan');
            data(j).I_Av = mean(intensities,'omitnan');
            data(j).I_Std = std(intensities,'omitnan');
            data(j).R2_Av = mean(R2,'omitnan');
            data(j).R2_Std = std(R2,'omitnan');
        %Else the average for the main ROI will be provided
        else
            data(j).CD_Av = particles;
            data(j).CD_Std = std(particles,'omitnan');
            data(j).BeamArea_Av = beamradius;
            data(j).BeamArea_Std = std(beamarea,'omitnan');
            data(j).w_Av = beamradius;
            data(j).w_Std = std(beamradius,'omitnan');
            data(j).I_Av = mean(intensities,'omitnan');
            data(j).I_Std = std(intensities,'omitnan');
            data(j).R2_Av = R2;
            data(j).R2_Std = std(R2,'omitnan');
        end
        
    %This is the non-scanning section   
    else
        %Generates a correlation function
        ICS2DCorr = corrfunc(croppedImageSeries);
        %Crops the correlation function to the region of interest
        ICS2DCorrCrop = autocrop(ICS2DCorr,12);
        %Fit of a 2D gaussian to the correlation function 
        %Parameters given in variable 'a'
        a = gaussfitold(ICS2DCorrCrop,'2d',pixelsize,wn,fitRadius);
        if strcmp(graphs,'y')
            figure(1)
            s=surf(ICS2DCorr(:,:,1));
            axis tight
            colormap(jet)
            xlabel('\eta','FontSize',12)
            ylabel('\xi','FontSize',12)
            zlabel('r(\xi,\eta)','FontSize',12)
            set(s,'LineStyle','none')
            title('Spatial Autocorrelation Function for First Image')
            
            figure(2)
            s=surf(ICS2DCorrCrop(:,:,1));
            axis tight
            colormap(jet)
            xlabel('\eta','FontSize',12)
            ylabel('\xi','FontSize',12)
            zlabel('r(\xi,\eta)','FontSize',12)
            set(s,'LineStyle','none')
            title('Cropped Spatial Autocorrelation Function for First Image')
            
            plotgaussfit(a(1,:),ICS2DCorrCrop(:,:,1),pixelsize,wn)
            hold off
        end
        %Generation of variables of interest from Gaussian fit
        particles = 1./a(:,1);
        beamArea = pi.*a(:,2).*a(:,3);
        beamRadius = (a(:,2)+a(:,3))./2;
        CD = particles./beamArea;
        
        %For no filtering selected, or only one gaussian was fit
        if any(1 == indx) | any(0 == indx) | length(CD) == 1
        
        %Filter options: Negative values -> 30% around w -> 3*std + mean
        else
            beamArea1 = beamArea;
            beamRadius1 = beamRadius;
            CD1 = CD;
            if ismember(4,indx)
                for i = 1:length(CD)
                    if beamArea1(i) < 0;
                        beamArea1(i) = NaN;
                    end
                    if beamRadius1(i) < 0;
                        beamRadius1(i) = NaN;
                    end
                    if CD1(i) < 0;
                        CD1(i) = NaN;
                    end
                end
            end
             if ismember(3,indx)
                for i = 1:length(CD)
                    if beamRadius1(i) < 0.7*realW | beamRadius(i) > 0.7*realW
                        beamRadius1(i) = NaN;
                        CD1(i) = NaN;
                        beamArea1(i) = NaN;
                    end
                end
            end
            if ismember(2, indx)
                baM = mean(beamArea1,'omitnan');
                brM = mean(beamRadius1,'omitnan');
                cdM = mean(CD1,'omitnan');
                baS = std(beamArea1,'omitnan');
                brS = std(beamRadius1,'omitnan');
                cdS = std(CD1,'omitnan');
                for i = 1:length(CD)
                    if beamArea1(i) > 3*baS + baM | beamArea1(i) < -3*baS + baM
                        beamArea1(i) = NaN;
                    end
                    if beamRadius1(i) > 3*brS + brM | beamRadius1(i) < -3*brS + brM
                        beamRadius1(i) = NaN;
                    end
                    if CD1(i) > 3*cdS + cdM | CD1(i) < -3*cdS + cdM
                        CD1(i) = NaN;
                    end
                end
            end
        end
        %Will provide you with all the extracted data from an ROI
        if strcmp(allV,'y')
            data(j).w_Av =  beamRadius;
            data(j).w_Std = std(beamRadius,'omitnan')'
            data(j).BeamArea_Av =  beamArea;
            data(j).BeamArea_Std = std(beamArea,'omitnan');
            data(j).CD_Av =  CD;
            data(j).CD_Std = std(CD,'omitnan');
        %Else will average the data from the ROI to provide a single result
        else
            data(j).w_Av = mean(beamRadius,'omitnan');
            data(j).w_Std = std(beamRadius,'omitnan')'
            data(j).BeamArea_Av = mean(beamArea,'omitnan');
            data(j).BeamArea_Std = std(beamArea,'omitnan');
            data(j).CD_Av = mean(CD,'omitnan');
            data(j).CD_Std = std(CD,'omitnan');
        end
        
        data(j).I_Av = mean2(croppedImageSeries);
        data(j).I_Std = std2(croppedImageSeries);
        data(j).X_Size = crop(1,3);
        data(j).Y_Size = crop(1,4);
    end
end