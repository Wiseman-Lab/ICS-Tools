function opt = setupSTICS(varargin)
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
    
    opt(n_runs) = struct();
    data = cell(n_runs, 2);
    
    for run = 1:n_runs
%% Load files        
        disp(['Setup for run #', num2str(run)]);
        opt(run).STICCS = input('Do cross-correlation? ');
        sameData = 0;
        if run > 1 && opt(run).STICCS == opt(run-1).STICCS
            sameData = input(['Use the same data as in run #', num2str(run-1), '? ']);
            if sameData
                opt(run).fileName{1} = opt(run-1).fileName{1}; 
                opt(run).fileName{2} = opt(run-1).fileName{2};
                opt(run).filePath{1} = opt(run-1).filePath{1};
                opt(run).filePath{2} = opt(run-1).filePath{2};
                opt(run).pixelSize =  opt(run-1).pixelSize;
                opt(run).timeFrame = opt(run-1).timeFrame;
                opt(run).CropArea = opt(run-1).CropArea;
                opt(run).maskCell = opt(run-1).maskCell;
            end
        end
        
        if ~isfield(opt(run), 'fileName') || isempty(opt(run).fileName) || ...
                ~isfield(opt(run), 'filePath') || isempty(opt(run).filePath)
            if opt(run).STICCS
                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 1 for run #', num2str(run)]);
                opt(run).fileName{1} = file_id(1:end-4);
                opt(run).filePath{1} = [file_path, file_id];
                stack1 = readFileToStack(opt(run).filePath{1});
                if size(stack1,1) ~= size(stack1,2)
                    stack1 = makeSquare(stack1);
                end
                data{run, 1}  = stack1;

                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 2 for run #', num2str(run)]);
                opt(run).fileName{2} = file_id(1:end-4);
                opt(run).filePath{2} = [file_path, file_id];
                stack2 = readFileToStack(opt(run).filePath{2});
                if size(stack2,1) ~= size(stack2,2)
                    stack2 = makeSquare(stack2);
                end                

                while size(stack1, 1) ~= size(stack2, 1) || size(stack1, 2) ~= size(stack2, 2)
                    disp('Error: Both series must be the same size!')
                    [file_id, file_path] = uigetfile("*", ['Select the file for Stack 2 for run #', num2str(run)]);
                    opt(run).fileName{2} = file_id(1:end-4);
                    opt(run).filePath{2} = [file_path, file_id];
                    stack2 = readFileToStack(opt(run).filePath{2});
                    if size(stack2,1) ~= size(stack2,2)
                        stack2 = makeSquare(stack2);
                    end                
                end
                data{run, 2} = stack2;

            else
                [file_id, file_path] = uigetfile("*", ['Select the file for Stack 1 for run #', num2str(run)]);
                opt(run).fileName{1} = file_id(1:end-4);
                opt(run).fileName{2} = opt(run).fileName{1};
                opt(run).filePath{1} = [file_path, file_id];
                opt(run).filePath{2} = opt(run).filePath{1};
                stack1 = readFileToStack(opt(run).filePath{1});
                if size(stack1,1) ~= size(stack1,2)
                    stack1 = makeSquare(stack1);
                end                
                stack2 = stack1;
                data{run,1} = stack1; data{run, 2} = stack2;
            end
            assignin("base", "data", data);            
        end
        
%% Metadata
% Manually input metadata. Upcoming options:
% Import metadata dfrom metadata.mat file
% Get metadata using getMetadata.m

        if ~sameData && ~exist('metadata', 'var')
            opt(run).pixelSize = input('Pixel size (um) = ');
            opt(run).timeFrame = input(['Interval between frames = ']);
            timeInSec = input('Is that in seconds? Otherwise, enter 0. ');
            if timeInSec
                opt(run).TimeUnits = 'sec';                
            else
                opt(run).TimeUnits = 'min';                                
            end
        elseif exist('metadata', 'var')
            
        end
        if opt(run).STICCS, opt(run).timeDelay = input('Time between channels (s)'); end

%% Crop image
        cropImage = input('Crop image? ');
        if cropImage %|| ~isfield(options(run), 'CropArea') || isempty(options(run).CropArea) 
            opt(run).OriginalSize(1,:) = size(stack1);
            opt(run).OriginalSize(2,:) = size(stack2);
            [stack1, cropArea] = serimcropold(stack1);
            opt(run).CropArea = cropArea;
            % Crop the same rectangle from channel 2
            for i = 1:size(stack2, 3)
                tempstack(:,:,i) = imcrop(stack2(:,:,i), cropArea); %#ok<AGROW>
                imwrite(stack1(:,:,i), [opt(run).filePath{1,1}(1,1:end-4), '_cropped4STICS.tif'], 'Compression', 'none', 'WriteMode', 'append');
                if opt(run).STICCS
                    imwrite(stack2(:,:,i), [opt(run).filePath{1,2}(1,1:end-4), '_cropped4STICS.tif'], 'Compression', 'none', 'WriteMode', 'append');
                end
            end
            data{run,1} = stack1; 
            data{run,2} = tempstack; clear tempstack;
            assignin("base", "data", data);            
        else
            opt(run).CropArea = [];
        end

%% Polygon mask        
        polygonMask = input('Do you want to select a polygon area? ');
        if polygonMask
            dynamicMask = input('Do you want to use a dynamic mask?');
            if dynamicMask
                [mask_Name, mask_Path] = uigetfile("*.mat", ['Select Mask File for run #', num2str(run), ', Channel 1']);
                load([mask_Path, mask_Name], 'maskFinal');
                
                % We're forcing the input stack to be square so we do the
                % same with the mask:
                if size(maskFinal{1,1},1) ~= size(maskFinal{1,1},2)
                    for i =1:size(maskFinal,2)
                        maskFinal{1,i} = makeSquare(maskFinal{1,i});
                        
                    end
                end
                
                % What is the size of the mask wrt the original and cropped
                % stack?
                [m1, m2] = size(maskFinal{1,1}); 
                
                o1 = opt(run).OriginalSize(1,1); 
                o2 = opt(run).OriginalSize(1,2); 

                stackSize = size(stack1); 
                s1 = stackSize(1, 1);
                s2 = stackSize(1, 2);
                s3 = stackSize(1, 3);
                
                if ~isequal([m1, m2], [o1, o2])
                    disp('Warning: Mask and Original Stack are not the same size.');
                    
                elseif isequal([m1, m2], [o1, o2])
                        % Mask and stack were the same size before
                        % cropping, so we can use the cropArea used to crop
                        % the stack if needed.                                                                            
                else
%                     marquee = 0;
%                     cropMask = 0;
                end
                
                if m1 > s1 && m2 > s2
                    % Mask is bigger than the cropped stack so 
                    % the mask needs to be cropped down to fit 
                    % the size of the stack.
                    cropMask =1; marquee = 0;
                    disp('Mask is bigger than stack');
                elseif s1 > m1 && s2 > m2
                    % The cropped stack is bigger than the mask so
                    % we need to center the mask wrt to the cropped
                    % stack based on the original size of the stack
                    cropMask = 0; marquee = 1;
                    disp('Stack is bigger than mask');
                end                 
                
                for  i = 1:s3
                    mask = maskFinal{1,i};
                    

                    
                    if marquee % Stack bigger than mask
                        if isempty(opt(run).CropArea)
                            % The stack has not been cropped so we can
                            % align the mask to the stack directly.
                            dynamicMask = zeros(size(stack1,1,2));
                            marginX = floor(s2 - m2)/2;
                            marginY = floor(s1 - m1)/2;
                            dynamicMask(marginY+1:end-marginY,...
                                marginX+1:end-marginX) = mask;
                        end
                    elseif isequal([m1, m2], [o1, o2]) && cropMask
                            % The stack has been cropped so we need to use
                            % its original size to center the mask before
                            % cropping it.
                            dynamicMask = zeros(o1, o2);
                            
                            marginX = floor((o2 - m2)/2);
                            marginY = floor((o1 - m1)/2);
                            
                            dynamicMask(marginY+1:end-marginY,...
                                marginX+1:end-marginX) = mask;
                            
                            croppedMask = imcrop(dynamicMask, opt(run).CropArea);
                            if sum(mask,'all') ~= sum(croppedMask,'all')
                                disp('Warning: The mask may have true values outside of the cropped area!');
                            else
                                dynamicMask = croppedMask;
                            end
                    else
                        dynamicMask = mask;
                    end                    
                    opt(run).dynamicMask(:,:,i) = dynamicMask;
                    imwrite(opt.dynamicMask(:,:,i), ...
                        [opt(run).filePath{1,1}(1,1:end-4), '_DynamicMask.tif'], ...
                        'Compression', 'none', 'WriteMode', 'append');
                end                
                opt(run).maskCell = logical(sum(opt(run).dynamicMask,3));                
            else % Polygon mask
                if opt(run).STICCS
                    opt(run).maskCell = selectCellMaskFromAverageImage(stack1, stack2);
                else
                    opt(run).maskCell = selectCellMaskFromAverageImage(stack1);
                end
                opt(run).dynamicMask = repmat(opt.maskCell, [1, 1, size(stack1, 3)]);
            end
        else
            opt(run).maskCell = ones(size(stack1, 1), size(stack1, 2));
            opt(run).dynamicMask = repmat(opt.maskCell, [1, 1, size(stack1,3)]);
        end       

%% Immobile fraction removal        
        immobileFilters = categorical({'FourierWhole','MovingAverage','butterIIR','none'});
        useFilter = input('Apply Immobile Filter?');
        if useFilter
            opt(run).filtering = input("Select an immobile substraction algorithm: 'FourierWhole','MovingAverage' or 'butterIIR': ");
        else 
            opt(run).filtering = 'none';
        end
        if strcmpi(opt(run).filtering, 'MovingAverage')
            opt(run).MoveAverage = input("Moving average window size (frames): ");
            while mod(opt(run).MoveAverage, 2) == 0
                disp('The moving average window size must be an odd number!');
                opt(run).MoveAverage = input("Moving average window size (frames): ");
            end
        end
        
        while ~iscategory(immobileFilters, opt(run).filtering)
            disp('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
            opt(run).filtering = input("Select an immobile substraction algorithm: 'FourierWhole','MovingAverage','butterIIR' or 'none': ");            
        end

%% Set parameters
% Displacement parameters
        opt(run).ROIsize = input(['ROI size for run #', num2str(run), ' = ']);
        opt(run).ROIshift = input(['ROI shift for run #', num2str(run), ' = ']);
        opt(run).TOIsize = input(['TOI size for run #', num2str(run), ' = ']);
        opt(run).TOIshift = input(['TOI shift for run #', num2str(run), ' = ']);        
        
% Fit options
        opt(run).tauLimit = input(['Tau limit (max = ' num2str(opt(run).TOIsize) ') for run #', num2str(run), ' = ']);% was 30
        opt(run).fitRadius = input(['Maximum fit radius (max = ', num2str(opt(run).ROIsize/2), ', units = px) for run #', num2str(run),' = ']);% was 32

% Filtering options
        opt(run).maxHalfWidth = input(['Max Width (suggestion = ', num2str(opt(run).pixelSize*opt(run).ROIsize/4), ') for run #', num2str(run), ' = ']);% was 100
        opt(run).threshVector = 5; %input(['STD vector threshold for run #', num2str(run), ' = ']);%was 16 Maximum orientation deviation (in standard deviations) among neighbouring vectors
        opt(run).maxV = Inf; %input(['Maximum Velocity permitted for run #', num2str(run), ' = ']);

% Plotting options
        titleSuggestion = [opt(run).fileName{1}, ...
            '_',num2str(opt(run).ROIsize), 'x', num2str(opt(run).ROIshift), '_', ...
            num2str(opt(run).TOIsize), 'x', num2str(opt(run).TOIshift), '_t', num2str(opt(run).tauLimit),  ...
            'r', num2str(opt(run).fitRadius), 'w', num2str(opt(run).maxHalfWidth ), ...
            'sd', num2str(opt(run).threshVector), 'v', num2str(opt(run).maxV)];
        acceptTitle = input(['Title: ', titleSuggestion, ' OK?']);
        if acceptTitle
            opt(run).axisTitle = titleSuggestion; display(['Axis title = ', titleSuggestion]);
            opt(run).outputName = titleSuggestion; display(['Output filename = ', titleSuggestion]);
        else
            opt(run).axisTitle = input('Axis title = ');
            opt(run).outputName = input('Output filename = ');            
        end
        opt(run).path = input('Output path = '); %
        
        TOImeanBackground = input('Use TOI mean for background?');
        if TOImeanBackground
            opt(run).bgImage = 'TOI mean';
        else
            disp('Please select a background image (''Original'', ''TOI mean'' or "Other")');
            opt(run).bgImage = input('Background image = ');
        end
        if strcmp(opt(run).bgImage, 'Other')
            [bg_file_id, bg_file_path] = uigetfile("*", ['Select the file to display in the background for run #', num2str(run)]);
            opt(run).bg_fileName{1} = bg_file_id(1:end-4);
            opt(run).bg_filePath{1} = [bg_file_path, bg_file_id];
            bg_series = readFileToStack(opt(run).bg_filePath{1});
            while size(stack1, 1) ~= size(bg_series, 1) || size(stack1, 2) ~= size(bg_series, 2) || size(stack1, 3) ~= size(bg_series, 3)
                disp('Error: Both series must be the same size!')
                [bg_file_id, bg_file_path] = uigetfile("*", ['Select the file to display in the background for run #', num2str(run)]);
                opt(run).bg_fileName{1} = bg_file_id(1:end-4);
                opt(run).bg_filePath{1} = [bg_file_path, bg_file_id];
                bg_series = readFileToStack(opt(run).bg_filePath{1});
            end
            assignin("base", "bg_series", bg_series);
        end

%% Dilate Cell Mask
        % Create a kernel or structuring element with a circle of ones that
        % has a radius equal to the ROI size so that at least one vector
        % more is calculated within the newly added area. Then dilate the 
        % cell mask with the structuring element. This step is added
        % to provide extra vectors in the periphery of the cell mask
        % to improve the accuracy of the interpolation.
        se = strel('disk', opt(run).ROIshift);
        opt(run).maskCell = imdilate(opt(run).maskCell, se);
        
   %% Define defaults in case they are missing (do not edit): Use isempty
    % instead of isfield? Create all fields in advance?
        if ~isfield(opt(run), 'pixelSize'), opt(run).pixelSize = 0.10 ; end % in micrometers
        if ~isfield(opt(run), 'timeFrame'), opt(run).timeFrame = 10; end %time delay (s) between subsequent frames
        if ~isfield(opt(run), 'timeDelay'), opt(run).timeDelay = 0; end %extra time delay (s) between two channels
        if ~isfield(opt(run), 'tauLimit'), opt(run).tauLimit = 20; end % what is highest time lag to compute stics of TOI
        %filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
        if ~isfield(opt(run), 'filtering'), opt(run).filtering ='FourierWhole'; end 
        if ~isfield(opt(run), 'MoveAverage'), opt(run).MoveAverage = 21; end
        if ~isfield(opt(run), 'fitRadius'), opt(run).fitRadius = 5; end  % how many pixels (radius) are considered around the peak of corr fn when fitting 2D Gaussian 
        if ~isfield(opt(run), 'maxHalfWidth'), opt(run).maxHalfWidth = 15; end % how big (pixels) do you allow the radius of corr fn to be...removes wide corr fns...
        if ~isfield(opt(run), 'threshVector'), opt(run).threshVector = 8; end % what is the threshold delta(v) (um/s) used in discarding spurious vectors
        %if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
        if ~isfield(opt(run), 'ROIsize'), opt(run).ROIsize = 16; end % ROI size in pixels
        if ~isfield(opt(run), 'ROIshift'), opt(run).ROIshift = 4; end %what is ROI centers shift in pixels
        if ~isfield(opt(run), 'TOIsize'), opt(run).TOIsize = 60; end % what is TOI size in frames
        if ~isfield(opt(run), 'TOIshift'), opt(run).TOIshift = 1; end %how much is TOI shifted
        if ~isfield(opt(run), 'axisTitle'), opt(run).axisTitle = [opt(run).fileName{1}, ...
            '_',num2str(opt(run).ROIsize), 'x', num2str(opt(run).ROIshift), '_', ...
            num2str(opt(run).TOIsize), 'x', num2str(opt(run).TOIshift), '_t', num2str(opt(run).tauLimit),  ...
            'r', num2str(opt(run).fitRadius), 'w', num2str(opt(run).maxHalfWidth ), ...
            'sd', num2str(opt(run).threshVector), 'v', num2str(opt.maxV)]; 
        end %what will be used on vector map title
        if ~isfield(opt(run), 'outputName'), opt(run).outputName = opt(run).axisTitle; end%'Velocity Map'; end %output file name
%         if ~isfield(options(run), 'path'), options(run).path ='C:\Users\YOGA-PW\Dropbox (Wiseman Research)\Data_Elvis\STICCSOutput'; end
        if ~isfield(opt(run), 'exportimages'), opt(run).exportimages ='n'; end %'y' if you want export a pdf of every vector map, 'n' if you do not
        if ~isfield(opt(run), 'imagesformat'), opt(run).imagesformat ='png'; end % 'png' or 'pdf' images...can add other option (formats) in plotSingleVectorMap code
        if ~isfield(opt(run), 'movieformat'), opt(run).movieformat ='mp4'; end % movie format can be avi, jpeg or mp4
        
%         assignin("base", "opt", options);
    end

elseif nargin ~= 0 && isstruct(varargin{1})
    opt = varargin{1};
%     make sure it has all the required fields though!!
%     requiredFields = 0; % A list of the required Fields
%     if fieldnames(options) ~= requiredFields
%         error('Required fields missing!')
%     end
end

end
        