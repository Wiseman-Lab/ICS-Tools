%Scans cropped image and performs multiple sICS measurements on the data
function [ICS2DCorr] = scanICS(croppedImageSeries, rect, fitRadius, pixelsize, ROIsize, realW, whitenoise, graphs, index, allV)
series = croppedImageSeries;
ICS2DCorr = struct('location_x',{},'location_y',{},'raw_data',{},'icorrfn',{},'autocrop',{}, 'gaussfit',{},'Particles',{},'BeamArea',{},'BeamRadius',{},'CD',{},'R2',{});
ROIshift = ROIsize/2;

% shift of ROI's amount & define the position of ROI
fracROIshift=ROIsize/ROIshift;
along_y = floor((size(series,2)-ROIsize)/ROIshift)+1;
along_x = floor((size(series,1)-ROIsize)/ROIshift)+1;
%Gives centre of each ROI shift
position_x  =   zeros(1,along_x);
position_y  =   zeros(1,along_y);
for i=1:along_x
    position_x(i) = 1/2 + (i-1)/fracROIshift*ROIsize + ROIsize/2;
end
position_x=repmat(position_x,along_y,1);
position_x=position_x';
for j=1:along_y
    position_y(j) = 1/2 + (j-1)/fracROIshift*ROIsize + ROIsize/2;
end
position_y=repmat(position_y,along_x,1);

%pre-process to determine if you want to analyse the whole FOV or select a
%polygon to analyze...in case of polygon, it will speed up the whole
%calculations because less ROI-TOI to consider...
% figure(1)
% imagesc(mean(series,3))
% prompt = {'want to select the polygon or use whole FOV?'};
%     dlg_title = 'y or n';
%     num_lines = 1;
%     def = {'y'};
%     answer = inputdlg(prompt,dlg_title,num_lines,def);
answer = 'n';
if strcmp(answer,'y')
    %MaskCell puts a polygon mask around the selected area
    %cropping the polygon you selected
    maskCell=[];
    while isempty(maskCell)
        [x,y,maskCell,poly1.x,poly1.y]=roipoly;
    end
elseif strcmp(answer,'n')
    maskCell=ones(size(series,1),size(series,2));
end
close
postouse=zeros(along_x,along_y);
for w=1:along_x
    for s=1:along_y
        indx=1+(w-1)/fracROIshift*ROIsize:(ROIsize*(1+((w-1)/fracROIshift)));
        indy=1+(s-1)/fracROIshift*ROIsize:(ROIsize*(1+((s-1)/fracROIshift)));
        indx = floor(indx);
        indy = floor(indy);
        if(maskCell(indx,indy)==1)
            postouse(w,s)=1;
        end
    end
end

[iroi jroi]=find(postouse);

for u=1:length(iroi) %accessing the data of only
    %the ROI that are defined within the selected polygon
    i=iroi(u);
    j=jroi(u);
    ICS2DCorr(u).location_x = rect(1) + position_x(i,1);
    ICS2DCorr(u).location_y = rect(2) + position_y(1,j);
    fprintf('analyzing ROI %i of %i. Date and time is %s \n',u, length(iroi),datestr(now));
    %define region to do ics on
    regionanalyze = series((1+(i-1)/fracROIshift*ROIsize):(ROIsize*(1+((i-1)/fracROIshift))), (1+(j-1)/fracROIshift*ROIsize):(ROIsize*(1+((j-1)/fracROIshift))),:);
    ICS2DCorr(u).raw_data = regionanalyze;
    icorrfn = corrfunc(regionanalyze);
    if strcmp(graphs,'y')
        figure(1)
        s=surf(icorrfn(:,:,1));
        axis tight
        colormap(jet)
        xlabel('\eta','FontSize',12)
        ylabel('\xi','FontSize',12)
        zlabel('r(\xi,\eta)','FontSize',12)
        set(s,'LineStyle','none')
        title('Spatial Autocorrelation Function for First Image')
    end
    ICS2DCorr(u).icorrfn = icorrfn;
    ICS2DCorr(u).autocrop = autocrop(ICS2DCorr(u).icorrfn, fitRadius);
    if strcmp(graphs,'y')
        figure(2)
        s=surf(ICS2DCorr(u).autocrop(:,:,1));
        axis tight
        colormap(jet)
        xlabel('\eta','FontSize',12)
        ylabel('\xi','FontSize',12)
        zlabel('r(\xi,\eta)','FontSize',12)
        set(s,'LineStyle','none')
        title('Cropped Spatial Autocorrelation Function for First Image')
    end
    
    [ICS2DCorr(u).gaussfit, ICS2DCorr(u).R2] = gaussfitold(ICS2DCorr(u).autocrop,'2d',pixelsize,whitenoise,fitRadius);
    if strcmp(graphs,'y')
        plotgaussfit(ICS2DCorr(u).gaussfit(1,:),ICS2DCorr(u).autocrop(:,:,1),pixelsize,whitenoise)
    end
end


b = zeros(size(ICS2DCorr,2),6,size(regionanalyze,3));
for i = 1:size(ICS2DCorr,2)
    for j = 1:6
        b(i,j,:) = ICS2DCorr(i).gaussfit(:,j);
    end
    if strcmp(allV,'n')
        ICS2DCorr(i).Particles = mean(1./b(i,1,:));
        ICS2DCorr(i).BeamArea = mean(pi.*b(i,2,:).*b(i,3,:));
        ICS2DCorr(i).BeamRadius =  mean(sqrt(b(i,2,:).*b(i,3,:)));
        ICS2DCorr(i).CD = mean((1./b(i,1,:))./(pi.*b(i,2,:).*b(i,3,:)));
    elseif strcmp(allV,'y')
        ICS2DCorr(i).Particles = (1./b(i,1,:));
        ICS2DCorr(i).BeamArea = (pi.*b(i,2,:).*b(i,3,:));
        ICS2DCorr(i).BeamRadius =  (sqrt(b(i,2,:).*b(i,3,:)));
        ICS2DCorr(i).CD = ((1./b(i,1,:))./(pi.*b(i,2,:).*b(i,3,:)));
    end
end

if any(1 == index) | any(0 == index)
else
    if ismember(4,index)
        for u = 1:size(ICS2DCorr,2)
            for k = 1:length(ICS2DCorr(u).CD)
                if ICS2DCorr(u).CD(k)<0
                    ICS2DCorr(u).CD(k) = NaN;
                end
                if ICS2DCorr(u).BeamArea(k)<0
                    ICS2DCorr(u).BeamArea(k) = NaN;
                end
                if ICS2DCorr(u).BeamRadius(k)<0
                    ICS2DCorr(u).BeamRadius(k) = NaN;
                end
            end
        end
    end
    if ismember(3,index)
        for u = 1:size(ICS2DCorr,2)
            for k = 1:length(ICS2DCorr(u).CD)
                if ICS2DCorr(u).BeamRadius(k) < 0.7*realW | ICS2DCorr(u).BeamRadius(k) > 1.3*realW
                    ICS2DCorr(u).CD(k) = NaN;
                    ICS2DCorr(u).BeamArea(k) = NaN;
                    ICS2DCorr(u).BeamRadius(k) = NaN;
                end
            end
        end
    end
    if ismember(2, index)
        values1 = 0;
        values2 = 0;
        values3 = 0;
        for u = 1:size(ICS2DCorr,2)
            values1 = [values1; ICS2DCorr(u).CD];
            values2 = [values2; ICS2DCorr(u).BeamArea];
            values3 = [values3; ICS2DCorr(u).BeamRadius];
        end
        values1 = values1(2:length(values1))
        values2 = values2(2:length(values2))
        values3 = values3(2:length(values3))
        for u = 1:size(ICS2DCorr,2)
            for k = 1:length(ICS2DCorr(u).CD)
                if ICS2DCorr(u).CD(k) > 3*std(values1,'omitnan')+mean(values1,'omitnan')| ICS2DCorr(u).CD(k) < -3*std(values1,'omitnan')+mean(values1,'omitnan')
                    ICS2DCorr(u).CD(k) = NaN;
                end
                if ICS2DCorr(u).BeamArea(k) > 3*std(values2,'omitnan')+mean(values2,'omitnan') |  ICS2DCorr(u).BeamArea(k) < -3*std(values2,'omitnan')+mean(values2,'omitnan')
                    ICS2DCorr(u).BeamArea(k) = NaN;
                end
                if ICS2DCorr(u).BeamRadius(k) > 3*std(values3,'omitnan')+mean(values3,'omitnan') |  ICS2DCorr(u).BeamRadius(k) < -3*std(values2,'omitnan')+mean(values2,'omitnan')
                    ICS2DCorr(u).BeamRadius(k) = NaN;
                end
            end
        end
    end
    
    
end



