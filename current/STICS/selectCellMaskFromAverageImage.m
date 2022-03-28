function [maskCell] = selectCellMaskFromAverageImage(series1, series2)
% function uses the average image of the series 
% and prompts user to select the
% polygon containing the cell. The ouputs, postouse and maskCell, can be
% saved and used later in STIC(C)S batch processing code.
% Inputs: path is the directory where files are loaded from
% name: name of the orignal file (string)
% series : 3D array (image series)
% ROIsize: size of the region of interest used in STICS, 16 for 16x16 ROI
% ROIshift: pixel amount by which ROI is shifted in x and y. For 16x16 set
% the default is 4...

figure;
if nargin == 2
    imcomp(mean(series1,3), mean(series2,3), 'y'); 
else
    imagesc(mean(series1,3)) 
end
axis image;

maskCell = [];
while isempty(maskCell)
    [~,~,maskCell,poly1.x,poly1.y] = roipoly;
end
close

end