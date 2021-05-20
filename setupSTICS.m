function options = setupSTICS(varargin)
%setupSTICS Creates an options structure with the parameters to be used for
%each STICS run. 

% if else loop to choose whether to use the GUI of the Command Window

%% Set batch mode
if nargin == 0
    % STICS options
    STICS_batch_mode = input('Do you want to analyze many files or parameter combinations in batch mode? ');
    if STICS_batch_mode
        n_runs = input('How many different runs (files / parameter combinations) do you want to execute? ');
    else
        n_runs = 1;
    end
    
    options(n_runs) = struct();
    data = cell(n_runs, 2);
    
    for run = 1:n_runs
%% Load files        
        disp(['Setup for run #', num2str(run)]);
        options(run).STICCS = input('Do cross-correlation? ');
        sameData = 0;
        if run > 1 && options(run).STICCS == options(run-1).STICCS
            sameData = input(['Use the same data as in run #', num2str(run-1), '? ']);
            if sameData
                options(run).fileName{1} = options(run-1).fileName{1}; 
                options(run).fileName{2} = options(run-1).fileName{2};
                options(run).filePath{1} = options(run-1).filePath{1};
                options(run).filePath{2} = options(run-1).filePath{2};
                options(run).pixelSize =  options(run-1).pixelSize;
                options(run).timeFrame = options(run-1).timeFrame;
                options(run).CropArea = options(run-1).CropArea;
                options(run).maskCell = options(run-1).maskCell;
            end
        end
        if ~isfield(options(run), 'fileName') || isempty(options(run).fileName) || ...
                ~isfield(options(run), 'filePath') || isempty(options(run).filePath)
            if options(run).STICCS
                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 1 for run #', num2str(run)]);
                options(run).fileName{1} = file_id(1:end-4);
                options(run).filePath{1} = [file_path, file_id];
                series1 = readFileToStack(options(run).filePath{1});
                data{run, 1}  = series1;
%                 assignin("base", "series1", series1);

                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 2 for run #', num2str(run)]);
                options(run).fileName{2} = file_id(1:end-4);
                options(run).filePath{2} = [file_path, file_id];
                series2 = readFileToStack(options(run).filePath{2});
%                 assignin("base", "series2", series2);

                while size(series1, 1) ~= size(series2, 1) || size(series1, 2) ~= size(series2, 2)
                    disp('Error: Both series must be the same size!')
                    [file_id, file_path] = uigetfile("*", ['Select the file for Stack 2 for run #', num2str(run)]);
                    options(run).fileName{2} = file_id(1:end-4);
                    options(run).filePath{2} = [file_path, file_id];
                    series2 = readFileToStack(options(run).filePath{2});
%                     assignin("base", "series2", series2);                
                end
                data{run, 2} = series2;
                assignin("base", "data", data);                

            else
                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 1 for run #', num2str(run)]);
                options(run).fileName{1} = file_id(1:end-4);
                options(run).fileName{2} = options(run).fileName{1};
                options(run).filePath{1} = [file_path, file_id];
                options(run).filePath{2} = options(run).filePath{1};
                series1 = readFileToStack(options(run).filePath{1});
                series2 = series1;
                data{run,1} = series1; data{run, 2} = series2;
                assignin("base", "data", data);
%                 assignin("base", "series1", series1);
            end
        end
        
%% Metadata
% Manually input metadata. Upcoming options:
% Import metadata dfrom metadata.mat file
% Get metadata using getMetadata.m

        if ~sameData
            options(run).pixelSize = input('Pixel size (um) = ');
            options(run).timeFrame = input('Interval between frames (s) = ');
        end
        if options(run).STICCS, options(run).timeDelay = input('Time between channels (s)'); end

%% Crop image
        cropImage = input('Crop image? ');
        if cropImage %|| ~isfield(options(run), 'CropArea') || isempty(options(run).CropArea) 
            [series1, cropArea] = serimcropold(series1);
%             assignin("base", "series1", series1);
            options(run).CropArea = cropArea;
            if options(run).STICCS
                % crop same rectangle within channel 2
                for i = 1:size(series2, 3)
                    series2(:,:,i) = imcrop(series2(:,:,i), cropArea);
                end
%                 assignin("base", "series2", series2);
            end
            data{run,1} = series1; data{run, 2} = series2;
            assignin("base", "data", data);            
        else
            options(run).CropArea = 0;
        end

%% Polygon mask        
        polygonMask = input('Do you want to select a polygon area? ');
        if polygonMask
            if options(run).STICCS
                options(run).maskCell = selectCellMaskFromAverageImage(series1, series2);
            else
                options(run).maskCell = selectCellMaskFromAverageImage(series1);
            end
        else
            options(run).maskCell = ones(size(series1, 1), size(series1, 2));
        end

%% Immobile fraction removal        
        immobileFilters = categorical({'FourierWhole','MovingAverage','butterIIR','none'});
        options(run).filtering = input("Select an immobile substraction algorithm: 'FourierWhole','MovingAverage','butterIIR' or 'none': ");
        if strcmpi(options(run).filtering, 'MovingAverage')
            options(run).MoveAverage = input("Moving average window size (frames): ");
            while mod(options(run).MoveAverage, 2) == 0
                disp('The moving average window size must be an odd number!');
                options(run).MoveAverage = input("Moving average window size (frames): ");
            end
        end
        
        while ~iscategory(immobileFilters, options(run).filtering)
            disp('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
            options(run).filtering = input("Select an immobile substraction algorithm: 'FourierWhole','MovingAverage','butterIIR' or 'none': ");            
        end

%% Set parameters
% Displacement parameters
        options(run).ROIsize = input(['ROI size for run #', num2str(run), ' = ']);
        options(run).ROIshift = input(['ROI shift for run #', num2str(run), ' = ']);
        options(run).TOIsize = input(['TOI size for run #', num2str(run), ' = ']);
        options(run).TOIshift = input(['TOI shift for run #', num2str(run), ' = ']);        
        
% Fit options
        options(run).tauLimit = input(['Tau limit for run #', num2str(run), ' = ']);% was 30
        options(run).fitRadius = input(['Maximum fit radius for run #', num2str(run),' = ']);% was 32

% Filtering options
        options(run).omegaThreshold = input(['Omega Threshold for run #', num2str(run), ' = ']);% was 100
        options(run).threshVector = input(['Orientation vector threshold for run #', num2str(run), ' = ']);%was 16 Maximum orientation deviation among neighbouring vectors

% Plotting options
        disp(['Suggestion: ', options(run).fileName{1}, ...
        '_',num2str(options(run).ROIsize), 'x', num2str(options(run).ROIshift), '_', ...
        num2str(options(run).TOIsize), 'x', num2str(options(run).TOIshift), '_t', num2str(options(run).tauLimit),  ...
        'r', num2str(options(run).fitRadius), 'o', num2str(options(run).omegaThreshold ), ...
        'v', num2str(options(run).threshVector)]);
        options(run).axisTitle = input('Axis title = ');
        options(run).outputName = input('Output filename = ');
        options(run).path = input('Output path = '); %
        disp('Please select a background image (''Original'', ''TOI mean'' or "Other")');
        options(run).bgImage = input('Background image = ');
        if strcmp(options(run).bgImage, 'Other')
            [bg_file_id, bg_file_path] = uigetfile("*", ['Select the file to display in the background for run #', num2str(run)]);
            options(run).bg_fileName{1} = bg_file_id(1:end-4);
            options(run).bg_filePath{1} = [bg_file_path, bg_file_id];
            bg_series = readFileToStack(options(run).bg_filePath{1});
            while size(series1, 1) ~= size(bg_series, 1) || size(series1, 2) ~= size(bg_series, 2) || size(series1, 3) ~= size(bg_series, 3)
                disp('Error: Both series must be the same size!')
                [bg_file_id, bg_file_path] = uigetfile("*", ['Select the file to display in the background for run #', num2str(run)]);
                options(run).bg_fileName{1} = bg_file_id(1:end-4);
                options(run).bg_filePath{1} = [bg_file_path, bg_file_id];
                bg_series = readFileToStack(options(run).bg_filePath{1});
            end
            assignin("base", "bg_series", bg_series);
        end


    % Define defaults in case they are missing (do not edit): Use isempty
    % instead of isfield? Create all fields in advance?
        if ~isfield(options(run), 'pixelSize'), options(run).pixelSize = 0.10 ; end % in micrometers
        if ~isfield(options(run), 'timeFrame'), options(run).timeFrame = 10; end %time delay (s) between subsequent frames
        if ~isfield(options(run), 'timeDelay'), options(run).timeDelay = 0; end %extra time delay (s) between two channels
        if ~isfield(options(run), 'tauLimit'), options(run).tauLimit = 20; end % what is highest time lag to compute stics of TOI
        %filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
        if ~isfield(options(run), 'filtering'), options(run).filtering ='FourierWhole'; end 
        if ~isfield(options(run), 'MoveAverage'), options(run).MoveAverage = 21; end
        if ~isfield(options(run), 'fitRadius'), options(run).fitRadius = 5; end  % how many pixels (radius) are considered around the peak of corr fn when fitting 2D Gaussian 
        if ~isfield(options(run), 'omegaThreshold'), options(run).omegaThreshold = 15; end % how big (pixels) do you allow the radius of corr fn to be...removes wide corr fns...
        if ~isfield(options(run), 'threshVector'), options(run).threshVector = 8; end % what is the threshold delta(v) (um/s) used in discarding spurious vectors
        %if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
        if ~isfield(options(run), 'ROIsize'), options(run).ROIsize = 16; end % ROI size in pixels
        if ~isfield(options(run), 'ROIshift'), options(run).ROIshift = 4; end %what is ROI centers shift in pixels
        if ~isfield(options(run), 'TOIsize'), options(run).TOIsize = 60; end % what is TOI size in frames
        if ~isfield(options(run), 'TOIshift'), options(run).TOIshift = 1; end %how much is TOI shifted
        if ~isfield(options(run), 'axisTitle'), options(run).axisTitle = [options(run).fileName{1}, ...
        '_',num2str(options(run).ROIsize), 'x', num2str(options(run).ROIshift), '_', ...
        num2str(options(run).TOIsize), 'x', num2str(options(run).TOIshift), '_t', num2str(options(run).tauLimit),  ...
        'r', num2str(options(run).fitRadius), 'o', num2str(options(run).omegaThreshold ), ...
        'v', num2str(options(run).threshVector)]; 
        end %what will be used on vector map title
        if ~isfield(options(run), 'outputName'), options(run).outputName = options(run).axisTitle; end%'Velocity Map'; end %output file name
%         if ~isfield(options(run), 'path'), options(run).path ='C:\Users\YOGA-PW\Dropbox (Wiseman Research)\Data_Elvis\STICCSOutput'; end
        if ~isfield(options(run), 'exportimages'), options(run).exportimages ='n'; end %'y' if you want export a pdf of every vector map, 'n' if you do not
        if ~isfield(options(run), 'imagesformat'), options(run).imagesformat ='png'; end % 'png' or 'pdf' images...can add other option (formats) in plotSingleVectorMap code
        if ~isfield(options(run), 'movieformat'), options(run).movieformat ='mp4'; end % movie format can be avi, jpeg or mp4
        
        assignin("base", "options", options);
    end

elseif nargin ~= 0 && isstruct(varargin{1})
    options = varargin{1};
%     make sure it has all the required fields though!!
%     requiredFields = 0; % A list of the required Fields
%     if fieldnames(options) ~= requiredFields
%         error('Required fields missing!')
%     end
end

end
        