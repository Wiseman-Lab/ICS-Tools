function [image_data] = readFileToStack(file_path)
%readFileToStack reads an image series of format TIF, LSM or MAT into a 3D
%stack.
%   Rodrigo Migueles Ramirez, March 2021.

% if data is in MultiTiff files
if strcmp(file_path(end-2:end),'tif')
    image_data = rd_img16(file_path);
% if data is in microscope formats such as Zeiss .lsm files, 'imreadBF' to
% load each channel data. Please open the 'imreadBF.m' and see how the
% input parameters are defined in order to load your data. 
elseif strcmp(file_path(end-2:end),'lsm')
    image_data = imreadBF(file_path, 1, 1:100, 1);
% if data is stored in a .mat variable
elseif strcmp(file_path(end-2:end),'mat')
    image_data=load(file_path);
else 
   error('Your input data set should be either ''.tif'' ,''.lsm'' or ''.mat'' '); 
end

% include ND2 and other formats using Bio-Formats?

end