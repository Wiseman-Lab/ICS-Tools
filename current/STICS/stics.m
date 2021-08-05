function [timecorr, varargout] = stics(varargin)

% July 10, 2003
% David Kolin
% Updated Nov 21, 2006 to allow for cross correlation
% Calculates the full time correlation function given 3D array of image series
% Usage:
% [timecorr] = stics(imgser,upperTauLimit) 
% OR
% [crossCorr] = stics(imgser1,imgser2,upperTauLimit)
% OR
% [crossCorr autoCorr1 autoCorr2] = stics(imgser1,imgser2,upperTauLimit)
% where timecorr is the spatio-temporal correlation function of imgser,
% calculated up to time lag upperTauLimit
% In the second case, timecorr is the spatio-temporal cross-correlation
% function of imgser1 and imgser2
% Modified April 13/07 to calculate xcorr 1 vs 2 and then 2 vs 1 and
% average them
% Updated July 15, 2009 to correct for edge effect

useWaitBar = 0;

if useWaitBar
    set(gcbf,'pointer','watch');
    h = waitbar(0,'Calculating time correlation functions...');
end

% tau is the lag
% pair is the nth pair of a lag time

if length(varargin)==2
    imgser = varargin{1};
    upperTauLimit = min(varargin{2},size(imgser,3));
    
    %timecorr = zeros(size(imgser,1),size(imgser,2),upperTauLimit);  % preallocates lagcorr matrix for storing raw time corr functions
    %SeriesMean = squeeze(mean(mean(imgser)));
    fftStack = double(zeros(size(imgser)));
    imageMean = squeeze(mean(mean(imgser,1),2));
    imageMean(find(imageMean==0))=10^-12;
    
    %{
    for i = 1:size(imgser,3)
        fftStack(:,:,i) = fft2(double(imgser(:,:,i)));
    end
    %}
    
    % Vectorized
    i=1:size(imgser,3);
    fftStack(:,:,i) = fft2(double(imgser(:,:,i)));
  
    for tau = 0:upperTauLimit-1
        lagcorr = zeros(size(imgser,1),size(imgser,2),(size(imgser,3)-tau));
        for pair=1:(size(imgser,3)-tau)
            %lagcorr(:,:,pair) = fftStack(:,:,pair).*conj(fftStack(:,:,(pair+tau)));
            lagcorr(:,:,pair) = ifft2(fftStack(:,:,pair).*conj(fftStack(:,:,(pair+tau))),'symmetric')./(imageMean(pair)*imageMean(pair+tau));
                
        end
        %timecorr(:,:,(tau+1)) = fftshift(ifft2mod(mean(lagcorr,3),'symmetric'));
        timecorr(:,:,(tau+1)) = fftshift(mean(lagcorr,3))./size(imgser,1)./size(imgser,2)-1;
                
        % Checks for significance of global maximum
         if (tau==(upperTauLimit-1))||(~correlationSignificance(timecorr(:,:,tau+1)))
             timecorr = timecorr(:,:,1:(end-1)); % cut off the "bad" lag
             break
         end
        
        if useWaitBar
            if ishandle(h)
                waitbar((tau+1)/(upperTauLimit),h)
            else
                break
            end
        end
    end

elseif length(varargin)==3
    imgser1 = varargin{1};
    imgser2 = varargin{2};
    upperTauLimit = min(varargin{3},size(imgser1,3));
    
    if size(imgser1)~=size(imgser2)
        error('Sizes of imgser1 and imgser2 must be equal.')
    end
    
    %timecorr = zeros(size(imgser1,1),size(imgser1,2),upperTauLimit);  % preallocates lagcorr matrix for storing raw time corr functions
    %SeriesMean = squeeze(mean(mean(imgser)));
    fftStack1 = double(zeros(size(imgser1)));
    fftStack2 = double(zeros(size(imgser2)));
    imageMean1 = squeeze(mean(mean(imgser1,1),2));
    imageMean2 = squeeze(mean(mean(imgser2,1),2));
    imageMean1(find(imageMean1==0))=10^-12; 
    imageMean2(find(imageMean2==0))=10^-12;
    % Calculates FFTs only once to save time
    % Might cause memory problems for huge data sets
    for i = 1:size(imgser1,3)
        fftStack1(:,:,i) = fft2(double(imgser1(:,:,i)));
        fftStack2(:,:,i) = fft2(double(imgser2(:,:,i)));
    end
    
    % Calculates cross correlation functions
    for tau = 0:upperTauLimit-1
        lagcorr = zeros(size(imgser1,1),size(imgser1,2),(size(imgser1,3)-tau));
        for pair=1:(size(imgser1,3)-tau)
            %lagcorr(:,:,pair) = fftStack1(:,:,pair).*conj(fftStack2(:,:,(pair+tau)));
            lagcorr(:,:,pair) = ifft2(fftStack1(:,:,pair).*conj(fftStack2(:,:,(pair+tau))),'symmetric')./(imageMean1(pair)*imageMean2(pair+tau));
            
        end
        %timecorr(:,:,(tau+1)) = fftshift(real(ifft2(mean(lagcorr,3))));
        timecorr(:,:,(tau+1)) = fftshift(mean(lagcorr,3))./size(imgser1,1)./size(imgser1,2)-1;
        
        % Checks for significance of global maximum
        if (tau==(upperTauLimit-1))||(~correlationSignificance(timecorr(:,:,tau+1)))
            timecorr = timecorr(:,:,1:(end-1)); % cut off the "bad" lag
            break
        end

        
        if useWaitBar
            if ishandle(h)
              waitbar((tau+1)/(upperTauLimit),h,'Calculating cross-correlation functions...')
            else
              break
            end
        end
    end
    
    % Does autocorrelations for each channel if there are appropriate
    % outputs
    if nargout == 3 || nargout==4
        
        
        
        %cross-corr 21
       for tau = 0:upperTauLimit-1
           lagcorr = zeros(size(imgser1,1),size(imgser1,2),(size(imgser1,3)-tau));
           for pair=1:(size(imgser1,3)-tau)
               %lagcorr(:,:,pair) = fftStack2(:,:,pair).*conj(fftStack1(:,:,(pair+tau)));
               lagcorr(:,:,pair) = ifft2(fftStack2(:,:,pair).*conj(fftStack1(:,:,(pair+tau))),'symmetric')./(imageMean2(pair)*imageMean1(pair+tau));
            
           end
           %timecorr21(:,:,(tau+1)) = fftshift(real(ifft2(mean(lagcorr,3))));
           timecorr21(:,:,(tau+1)) = fftshift(mean(lagcorr,3))./size(imgser1,1)./size(imgser1,2)-1;
        
           % Checks for significance of global maximum for cross-corr 21
           if (tau==(upperTauLimit-1))||(~correlationSignificance(timecorr21(:,:,tau+1)))
              timecorr21 = timecorr21(:,:,1:(end-1)); % cut off the "bad" lag
              break
           end
           if useWaitBar
                if ishandle(h)
                  waitbar((tau+1)/(upperTauLimit),h,'Calculating auto-correlation functions...')
                else
                  break
                end
           end
        end   
        
        
        
        for tau = 0:upperTauLimit-1
            lagcorr = zeros(size(imgser1,1),size(imgser1,2),(size(imgser1,3)-tau));
            for pair=1:(size(imgser1,3)-tau)
                %lagcorr(:,:,pair) = fftStack1(:,:,pair).*conj(fftStack1(:,:,(pair+tau)));
                lagcorr(:,:,pair) = ifft2(fftStack1(:,:,pair).*conj(fftStack1(:,:,(pair+tau))),'symmetric')./(imageMean1(pair)*imageMean1(pair+tau));
            
            end
            %timecorr1(:,:,(tau+1)) = fftshift(real(ifft2(mean(lagcorr,3))));
            timecorr1(:,:,(tau+1)) = fftshift(mean(lagcorr,3))./size(imgser1,1)./size(imgser1,2)-1;
        
            % Checks for significance of global maximum
            if (tau==(upperTauLimit-1))||(~correlationSignificance(timecorr1(:,:,tau+1)))
                timecorr1 = timecorr1(:,:,1:(end-1)); % cut off the "bad" lag
                break
            end

            if useWaitBar
                if ishandle(h)
                  waitbar((tau+1)/(upperTauLimit),h,'Calculating auto-correlation functions...')
                else
                  break
                end
            end
        end
        
        
        
        for tau = 0:upperTauLimit-1
            lagcorr = zeros(size(imgser1,1),size(imgser1,2),(size(imgser1,3)-tau));
            for pair=1:(size(imgser1,3)-tau)
                lagcorr(:,:,pair) = fftStack2(:,:,pair).*conj(fftStack2(:,:,(pair+tau)));
                lagcorr(:,:,pair) = ifft2(fftStack2(:,:,pair).*conj(fftStack2(:,:,(pair+tau))),'symmetric')./(imageMean2(pair)*imageMean2(pair+tau));
            end
            %timecorr2(:,:,(tau+1)) = fftshift(real(ifft2(mean(lagcorr,3))));
            timecorr2(:,:,(tau+1)) = fftshift(mean(lagcorr,3))./size(imgser2,1)./size(imgser2,2)-1;
        
            % Checks for significance of global maximum
            if (tau==(upperTauLimit-1))||(~correlationSignificance(timecorr2(:,:,tau+1)))
                timecorr2 = timecorr2(:,:,1:(end-1)); % cut off the "bad" lag
                break
            end

            if useWaitBar
                if ishandle(h)
                  waitbar((tau+1)/(upperTauLimit),h,'Calculating auto-correlation functions...')
                else
                  break
                end
            end
        end
        varargout(1) = {timecorr21};
        varargout(2) = {timecorr1};
        varargout(3) = {timecorr2};
    end
else
    error('Number of input arguments must be either 2 or 3.  See function help.')
end

if useWaitBar
    close(h)
end
set(gcbf,'pointer','arrow');