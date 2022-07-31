function [dataFiltered]= butterIIR(image_data, fc, extension,type)

% Time filtering using a highpass Butterworth IIR filter
% July 8th 2011
% Laurent Potvin-Trottier
% 
% image_data is the image time series 3D matrix
% fc is the normalized frequency cutoff (between 0 and 0.5)
%
% fc is normalized by the sampling rate so that the frequency cutoff in Hz
% is fc*f_sampling
%
% The function uses the most precise way to convert an analog 3rd order
% Butterworth IIR filter into a digital one and the apply it 2 times in
% forward and reverse order in each pixel of image_data. This create a zero
% phase 6th order filter.

% 
%fc=cutoffVelocity/60/pixelSize/6/(samplingRate/2); or something like
%that..
if nargin ==4
    filtertype=type;
else
    filtertype='high';
end

    
[z,p,k] = butter(1,fc,filtertype);
[sos,g] = zp2sos(z,p,k);      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);   % Create a dfilt object

% Plot the frequency response of the squared filter
%  hc = cascade(Hd,Hd);
%  figure;fvtool(hc);
if nargin ==3
    if strcmp(extension,'periodic')
        data=zeros(size(image_data).*[1 1 3]);
        image_data=double(image_data);
        % Periodic extension with adjusted means to allow for a small trend
        data(:,:,1:size(image_data,3))=image_data+repmat(mean(image_data(:,:,1:5),3)-mean(image_data(:,:,end-5:end),3),[1 1 size(image_data,3)]);
        data(:,:,1+size(image_data,3):2*size(image_data,3))=image_data;
        data(:,:,2*size(image_data,3)+1:end)=image_data-repmat(mean(image_data(:,:,1:5),3)-mean(image_data(:,:,end-5:end),3),[1 1 size(image_data,3)]);
      
        dataFiltered=filter(Hd,flipdim(filter(Hd,flipdim(data,3),3),3),3);
        dataFiltered=dataFiltered(:,:,1+size(image_data,3):2*size(image_data,3));
    else
        if strcmp(extension,'mirror')
            
        else 
            
            if strcmp(extension,'none')
                   dataFiltered=filter(Hd,flipdim(filter(Hd,flipdim(double(image_data),3),3),3),3);
         
            else
            error('Please select a valid extension (''periodic'' ,''mirror'' or ''none'')');
            end
        end
        
    end
else
    % Reverse-filter-reverse-filter the data
    dataFiltered=filter(Hd,flipdim(filter(Hd,flipdim(image_data,3),3),3),3);
end


