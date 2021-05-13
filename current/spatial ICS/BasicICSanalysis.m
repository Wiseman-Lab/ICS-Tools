%% Basic ICS, TICS and STICS
image_data = croppedImageSeries;
%% 
choice = 'Image Correlation Spectroscopy (ICS) for <N>'%chooseMethod;
if strcmp(choice,'Image Correlation Spectroscopy (ICS) for <N>')
    prompt = {'Pixel size (um)','Analyse image #','Fit Average (y or n)','Exclude G(0,0) from fit (y or n)?'};
dlg_title = 'Enter analysis Parameters';
defaultans = {num2str(0.1),num2str(1),'n','n'};
answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
pixelsize=str2double(answer{1});
im=str2double(answer{2});
analyseAve=answer{3};
removeG00=answer{4};
% calculate ICS 
ICS2DCorr = corrfunc(image_data);
ICS2DCorrCrop = autocrop(ICS2DCorr,20);
if strcmp(analyseAve,'n')
   a  = gaussfitold(ICS2DCorrCrop(:,:,im),'2d',pixelsize,removeG00);
   subplot(1,2,1)
   s=surf(ICS2DCorrCrop(:,:,im));
   axis tight
   colormap(jet)
   xlabel('\eta','FontSize',12)
   ylabel('\xi','FontSize',12)
   zlabel('r(\xi,\eta)','FontSize',12)
   set(s,'LineStyle','none')
   title(['Cropped Spatial Autocorrelation Function for Image' num2str(im)])
   set(gcf,'Color','black') 
   subplot(1,2,2)
   title(['Fit for Cropped Spatial Autocorrelation Function for Image' num2str(im)])
    plotgaussfit(a(1,1:6),ICS2DCorrCrop(:,:,im),pixelsize,removeG00);
    set(gcf,'Position',[300 300 1015 466])
    % results of fit
    particlesPerBeamArea = 1/a(1,1)
    beamArea = pi*a(:,2)*a(:,3)
    density = particlesPerBeamArea/beamArea
     zlim=get(gca,'ZLim');
    xlim=get(gca,'XLim');
    ylim=get(gca,'YLim');
    text(xlim(2),ylim(2),0.9*zlim(2),['Beam Area =' num2str(beamArea) ' \mum^{2}'],'Color',[1 1 1])
    text(xlim(2),ylim(2),0.8*zlim(2),['Particles Per Beam Area =' num2str(particlesPerBeamArea)],'Color',[1 1 1])
    text(xlim(2),ylim(2),0.7*zlim(2),['Particles Density =' num2str(density) ' \mum^{-2}'],'Color',[1 1 1])
    
elseif strcmp(analyseAve,'y')
    aveCropCorr=mean(ICS2DCorrCrop,3);
     a  = gaussfitold(aveCropCorr,'2d',pixelsize,removeG00);
   subplot(1,2,1)
   s=surf(aveCropCorr);
   axis tight
   colormap(jet)
   xlabel('\eta','FontSize',12)
   ylabel('\xi','FontSize',12)
   zlabel('r(\xi,\eta)','FontSize',12)
   set(s,'LineStyle','none')
   title('Average Cropped Spatial Autocorrelation Function','Color',[1 1 1])
   set(gcf,'Color','black') 
   subplot(1,2,2)
   title('Fit of Average Cropped Spatial Autocorrelation Function','Color',[1 1 1])
    plotgaussfit(a(1,1:6),aveCropCorr,pixelsize,removeG00);
    set(gcf,'Position',[300 300 1015 466])
    % results of fit
    particlesPerBeamArea = 1/a(1,1)
    beamArea = pi*a(:,2)*a(:,3)
    density = particlesPerBeamArea/beamArea
    zlim=get(gca,'ZLim');
    xlim=get(gca,'XLim');
    ylim=get(gca,'YLim');
    text(xlim(2),ylim(2),0.9*zlim(2),['Beam Area =' num2str(beamArea) ' \mum^{2}'],'Color',[1 1 1])
    text(xlim(2),ylim(2),0.8*zlim(2),['Particles Per Beam Area =' num2str(particlesPerBeamArea)],'Color',[1 1 1])
    text(xlim(2),ylim(2),0.7*zlim(2),['Particles Density =' num2str(density) ' \mum^{-2}'],'Color',[1 1 1])
    
else
     error('Answer should either be ''y'' or ''n''.')
end
elseif strcmp(choice,'Temporal-ICS for speed and diffusion') % TICS menu
    close
    GtDiff = tics(image_data,1);
    plot(GtDiff(:,1),GtDiff(:,2),'o')
    set(gca,'XScale','log');
    xlabel('\tau (frames)','FontSize',20)
    ylabel('<r_{11}(0,0,\tau)>','FontSize',20)
    set(gcf,'Color',[1 1 1])
    set(gca,'FontSize',20)
    %choicetics = menu('Choose which model to fit TICS fn ','1 component Diffusion','1 component flow','Sum of 1 component diffusing and 1 component flowing','Single component 3D diffusion');
    choicetics = chooseFitting;
    close 
     prompt = {'Frame Time (s)','PSF e^{-2} radius (um)','tau_{max} to fit (frames)','Axial PSF e^{-2} radius (um)'};
     dlg_title = 'Enter analysis Parameters';
     defaultans = {num2str(1),num2str(0.4),num2str(round(size(image_data,3)/5)),num2str(0)};
     answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
     timesize=str2double(answer{1});
     PSFSize=str2double(answer{2});
     taumax=str2double(answer{3});
     PSFZ=str2double(answer{4});
    if strcmp(choicetics,'1 component Diffusion')
    diffFitting = difffit(GtDiff(1:taumax,1)*timesize,GtDiff(1:taumax,2),PSFSize); 
     elseif strcmp(choicetics,'1 component flow')
    flowFitting = flowfit(GtDiff(1:taumax,1)*timesize,GtDiff(1:taumax,2),PSFSize);
     elseif strcmp(choicetics,'Sum of 1 component diffusing and 1 component flowing')
     diffflowFitting= diffflowfit(GtDiff(1:taumax,1)*timesize,GtDiff(1:taumax,2),PSFSize);   
     elseif strcmp(choicetics,'Single component 3D diffusion')
     Dfif3dFitting = difffit3d(GtDiff(1:taumax,1)*timesize,GtDiff(1:taumax,2),PSFSize,PSFZ);    
    end
elseif strcmp(choice,'Spatio-Temporal-ICS for flow')
    close
    prompt = {'Frame Time (s)','Pixel size (um)','tau_{max} to fit (frames)','Immobile filtering (y or n)'};
     dlg_title = 'Enter analysis Parameters';
     defaultans = {num2str(1),num2str(0.1),num2str(round(size(image_data,3)/5)),'n'};
     answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
     timesize=str2double(answer{1});
     pixelsize=str2double(answer{2});
     taumax=str2double(answer{3});
     filtering=answer{4};
     if strcmp(filtering,'y');
         image_data=immfilter(image_data);
     end
     corr=stics_byfft(image_data,mean(image_data(:)),taumax);
     sv(corr,'c',5)
    [Vx,Vy] = velocity(image_data,timesize,pixelsize,filtering,taumax,'n');
    set(gcf,'Position',[1043         350         560         420])
elseif strcmp(choice,'FCS')
  close
%  choice
    [corr,lags]=fcs_byfft(image_data,round(length(image_data)/100));
    plot(lags,corr,'o')
    set(gca,'XScale','log');
    xlabel('\tau (time points)','FontSize',20)
    ylabel('<G(\tau)>','FontSize',20)
    set(gcf,'Color',[1 1 1])
    set(gca,'FontSize',20)
    %choicetics = menu('Choose which model to fit FCS fn ','1 component Diffusion','1 component flow','Sum of 1 component diffusing and 1 component flowing','Single component 3D diffusion');
    choicefcs = chooseFitting;
    close 
     prompt = {'Sampling time (s)','PSF e^{-2} radius (um)','tau_{max} to fit','Axial PSF e^{-2} radius (um)'};
     dlg_title = 'Enter analysis Parameters';
     defaultans = {num2str(1),num2str(0.4),num2str(round(length(image_data)/100)),num2str(0)};
     answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
     timesize=str2double(answer{1});
     PSFSize=str2double(answer{2});
     taumax=str2double(answer{3});
     PSFZ=str2double(answer{4});
      if strcmp(choicefcs,'1 component Diffusion')
    diffFitting = difffit_fcs(lags(1:taumax)*timesize,corr(1:taumax),PSFSize,PSFZ); 
     elseif strcmp(choicefcs,'1 component flow')
    flowFitting = flowfit_fcs(lags(1:taumax)*timesize,corr(1:taumax),PSFSize,PSFZ);
     elseif strcmp(choicefcs,'Sum of 1 component diffusing and 1 component flowing')
     diffflowFitting= diffflowfit_fcs(lags(1:taumax)*timesize,corr(1:taumax),PSFSize,PSFZ);   
     elseif strcmp(choicefcs,'Single component 3D diffusion')
     Diff3dFitting = difffit3d_fcs(lags(1:taumax)*timesize,corr(1:taumax),PSFSize,PSFZ);    
    end
end