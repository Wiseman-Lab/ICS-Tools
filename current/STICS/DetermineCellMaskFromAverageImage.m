function [postouse,maskCell]=DetermineCellMaskFromAverageImage(path,name,series,ROIsize,ROIshift);

% function uses the average image of the series 
% and prompts user to select the
% polygon containing the cell. The ouputs, postouse and maskCell, can be
% saved and used later in STIC(C)S batch processing code.

%Inputs: path is the directory where files are loaded from
% name: name of the orignal file (string)
% series : 3D array (image series)
% ROIsize: size of the region of interest used in STICS, 16 for 16x16 ROI
% ROIshift: pixel amount by which ROI is shifted in x and y. For 16x6 set
% the default is 4...

along_y = floor((size(series,2)-ROIsize)/ROIshift)+1;
along_x = floor((size(series,1)-ROIsize)/ROIshift)+1;
fracROIshift=ROIsize/ROIshift;

for i=1:along_x
position_x(i) = 1/2 + (i-1)/fracROIshift*ROIsize + ROIsize/2;
end
position_x=repmat(position_x,along_y,1);
position_x=position_x';
for j=1:along_y
position_y(j) = 1/2 + (j-1)/fracROIshift*ROIsize + ROIsize/2;
end
position_y=repmat(position_y,along_x,1);

imagesc(mean(series,3)) 
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
    maskCell=ones(size(series,1),size(series,2));
end
close
postouse=zeros(size(position_x,1),size(position_x,2));
for w=1:size(position_x,1)
for s=1:size(position_x,2)
coorx=round(position_x(w,s,1));
coory=round(position_y(w,s,1));
if(maskCell(coorx,coory)==1)
postouse(w,s)=1;
end
end
end

save([path 'MaskWithinCell_Serie' name '.mat'],'postouse','maskCell')
end

