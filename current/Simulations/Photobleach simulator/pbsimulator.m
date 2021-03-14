function [imageseries, pixelsize, timesize] = pbsimulator()
%Kobi Pollard 2021

 prompt = {'Total Size (pixels)','Pixel Size (um)','Frames','Timesize (s)','Bleaching Constant (k)','Densities (arbitrary units)','Oligimerization (comma separated)','PSF (um)'};
    dlg_title = 'Enter Simulation Parameters';
    defaultans = {num2str(200),num2str(0.08125),num2str(100),num2str(0.05),num2str(0.05),num2str(0),num2str(0),num2str(0.3)};
    answer = inputdlg(prompt,dlg_title,[1 60],defaultans);
    totalsize = str2double(answer{1});
    pixelsize =str2double(answer{2});
    frames = str2double(answer{3});
    timesize = str2double(answer{4});
    k = str2double(answer{5});
    Density =str2num(answer{6});
    agg = str2num(answer{7});
    psf = str2double(answer{8});

%Density =~ den/radius^2 * 0.42
%Density = [40 40 40]; %/um^2
%agg = [1 2 3]; %-mer

radius = round(psf/pixelsize);
den = Density./radius^2 .* 0.42; 
extra = 1;
 
[matrix_final,particleimage_out,NUM,matrix_conv,matrix1] = create_and_convolve(totalsize,den,agg,radius,extra);

particles = sum(sum(particleimage_out));


[imageseries] = photobleach(matrix_conv,particleimage_out,matrix1, radius, timesize, k, frames, extra);


sv(imageseries,'c',5)

