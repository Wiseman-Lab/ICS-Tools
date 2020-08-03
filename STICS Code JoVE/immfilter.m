function movie = immfilter(movie,filterType,seriesNumber)

set(gcbf,'pointer','watch');
try
    h = timebar(['Immobile filtering image series ' num2str(seriesNumber) '...']);
catch % for backwards compatibility in case seriesNumber wasn't supplied
    h = timebar('Immobile filtering ...');
end


if strcmp(filterType,'Mobile');
    timebar(h,1/3)
    movie = fft(double(movie),[],3);
    timebar(h,2/3)
    avPixInt = mean(mean(movie(:,:,1)))/size(movie,3);
    movie(:,:,2:end) = 0;
    movie = uint16(real(ifft(movie,[],3))+avPixInt);
% 'F' option FFT filters the whole movie at once
% fastest for regular sized movies
elseif strcmp(filterType,'F')
	if ~ishandle(h); warndlg('Image series partially filtered.','Warning!'); set(gcbf,'pointer','arrow'); return; end
    timebar(h,1/3)
    movie = fft(double(movie),[],3);
    if ~ishandle(h); warndlg('Image series partially filtered.','Warning!'); set(gcbf,'pointer','arrow'); return; end
    timebar(h,2/3)
    avPixInt = mean(mean(movie(:,:,1)))/size(movie,3);
    movie(:,:,1) = 0;
    movie = uint16(real(ifft(movie,[],3))+avPixInt);
    % 'FM' option FFT filters each pixel location
    % one at a time.  Use for large movies... uses less
    % memory than 'F' option, but probably a little slower
elseif strcmp(filterType,'FM')
    AvIntSum = 0;
    for i = 1:size(movie,1)
        if ~ishandle(h); warndlg('Image series partially filtered.','Warning!'); set(gcbf,'pointer','arrow'); return; end
        timebar(h,i/size(movie,1))
        for j = 1:size(movie,2)
            fftPixelTrace = fft(double(movie(i,j,:)));
            AvIntSum = AvIntSum + fftPixelTrace(1);
            fftPixelTrace(1) = 0;
            movie(i,j,:) = int16(ifft(fftPixelTrace));
        end
    end
    movie = movie + AvIntSum/numel(movie);
    % If filtertype is a number, than use this as a moving average window size
elseif strcmp(filterType,'F0')
    for i = 1:size(movie,1)
        if ~ishandle(h); warndlg('Image series partially filtered.','Warning!'); set(gcbf,'pointer','arrow'); return; end
        timebar(h,i/size(movie,1))
        for j = 1:size(movie,2)
            fftPixelTrace = fft(double(movie(i,j,:)));
            fftPixelTrace(1) = 0;
            movie(i,j,:) = int16(ifft(fftPixelTrace));
        end
    end
    % If filtertype is a number, than use this as a moving average window
    % size
elseif isnumeric(filterType)
    % This old way is commented out
    % old way #1 -- using matlab's filter
    %     movieMean = filter(ones(1,filterType)/filterType,1,double(movie),[],3);
    %     movie = int16(double(movie) - movieMean);

    % old way number #2 -- using moving_average
    F = (filterType - 1)/2;
    for i = 1:size(movie,1)
        if ~ishandle(h); warndlg('Image series partially filtered.','Warning!'); set(gcbf,'pointer','arrow'); return; end
        timebar(h,i/size(movie,1))
        for j = 1:size(movie,2)
            averagePixelTrace = uint16(moving_average(double(movie(i,j,:)),F));
            movie(i,j,:) = movie(i,j,:)-averagePixelTrace;
        end
    end
elseif strcmp(filterType,'none')
    % do nothing
else
    error('Filtertype must be ''F'', ''FM'', ''none'', or for a moving average, a number.');
end

if ishandle(h)
    close(h)
end
set(gcbf,'pointer','arrow');
