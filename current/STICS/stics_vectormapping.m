function [velocityMap, opt] = stics_vectormapping(opt)
%this functions takes the timeseries and calculates for defined regions of
%size 'ROIsize' and 'TOIsize' stics correlation function and fits symmetric
%Gaussian to extract the flow in x and y directions...
%Developed by Elvis Pandzic 2012-2014 under supervision of 
%Prof. Paul W. Wiseman @ McGill University. 
% correspondence addressed to : e.pandzic@unsw.edu.au


    % Look at opt structure and either assign default value, or use passed
    % value
    if nargin==1
%         opt.pixelSize = 0.1; % just initialize it to something so the code below doesn't throw an error
        
    end
    
    if ~isfield(opt, 'pixelSize'), opt.pixelSize = 1; end %pixel size in um
    if ~isfield(opt, 'timeFrame'), opt.timeFrame = 1; end %time delay (s) between subsequent frames
    if ~isfield(opt, 'TimeUnits'), opt.TimeUnits = 'sec'; end
    if ~isfield(opt, 'VelUnits'), opt.VelUnits = 'um/sec'; end
    try
        stack1 = readFileToStack(opt.filePath{1,1});
        stack2 = readFileToStack(opt.filePath{1,2});
    catch

        [fileID, filePath] = uigetfile("*.tif", 'Could not find the requested file, please select it manually. ');
        stack1 = readFileToStack([filePath, fileID]);

        [fileID, filePath] = uigetfile("*.tif", 'Could not find the requested file, please select it manually. ');
        stack2 = readFileToStack([filePath, fileID]);
    end
        
    
    while size(stack1) ~= size(stack2)
        [fileID, filePath] = uigetfile("*.tif", 'Stacks must be the same size. Please select Stack 1 again. ');
        stack1 = readFileToStack([filePath, fileID]);
        
        [fileID, filePath] = uigetfile("*.tif", 'Stacks must be the same size. Please select Stack 2 again. ');
        stack2 = readFileToStack([filePath, fileID]);        
    end
    
    if ~isempty(opt.CropArea)
        cropped1 = zeros(size(opt.maskCell,1), size(opt.maskCell,2), size(stack1, 3));
        cropped2 = zeros(size(opt.maskCell,1), size(opt.maskCell,2), size(stack1, 3));
        for frame = 1:size(stack1,3)
            cropped1(:,:,frame) = imcrop(stack1(:,:,frame), opt.CropArea); 
            cropped2(:,:,frame) = imcrop(stack2(:,:,frame), opt.CropArea); 
        end
        clear stack1; stack1= cropped1; clear cropped1;
        clear stack2; stack2= cropped2; clear cropped2;
    end
    
    if ~isfield(opt, 'tauLimit'), opt.tauLimit = size(stack1,3); end % what is highest time lag to compute stics of TO
    %filtering: choose 'FourierWhole','MovingAverage','butterIIR','none'
    if ~isfield(opt, 'filtering'), opt.filtering = 'none'; end

    if ~isfield(opt, 'ROIsize'), opt.ROIsize = 32; end % ROI size in pixels
    if ~isfield(opt, 'ROIshift'), opt.ROIshift = 8; end %what is ROI centers shift in pixels
    if ~isfield(opt, 'TOIsize'), opt.TOIsize = 5; end % what is TOI size in frames
    if ~isfield(opt, 'TOIshift'), opt.TOIshift = 1; end %how much is TOI shifted
    
    if ~isfield(opt, 'maxHalfWidth'), opt.maxHalfWidth = ((opt.ROIsize/2)*opt.pixelSize)/2; end % % how big (um) do you allow the radius of corr fn to be...removes wide corr fns... % AKA omegaThreshold
    if ~isfield(opt, 'MoveAverage'), opt.MoveAverage = 21; end
    if ~isfield(opt, 'fitRadius'), opt.fitRadius = opt.ROIsize/4; end %how many pixels (radius) are considered around the peak of corr fn when fitting 2D Gaussian
    if ~isfield(opt, 'threshVector'), opt.threshVector = 5; end % what is the threshold delta(v) (um/s) used in discarding spurious vectors
    %if ~isfield(opt, 'threshVectorDotProd'), opt.threshVectorDotProd = 0.5; end
    if ~isfield(opt, 'maxV'); opt.maxV = Inf;end %((opt.ROIsize*opt.pixelSize/3)/opt.timeFrame)*60;end
    if ~isfield(opt, 'maskCell'), opt.maskCell = ones(size(stack1, 1), size(stack1, 2)); end

    if ~isfield(opt, 'path'), opt.path = cd; end
    if ~isfield(opt, 'axisTitle'), opt.axisTitle = [opt.fileName{1}, ...
        '_',num2str(opt.ROIsize), 'x', num2str(opt.ROIshift), '_', ...
        num2str(opt.TOIsize), 'x', num2str(opt.TOIshift), '_t', num2str(opt.tauLimit),  ...
        'r', num2str(opt.fitRadius), 'w', num2str(opt.maxHalfWidth ), ...
        'sd', num2str(opt.threshVector), 'v', num2str(opt.maxV)]; 
    end
    if ~isfield(opt, 'outputName'), opt.outputName = opt.axisTitle; end %output file name
    if ~isfield(opt, 'CheckSignificance'), opt.CheckSignificance = 1; end
    if ~isfield(opt, 'AllFrames'), opt.AllFrames = 1; end
    if ~isfield(opt, 'FrameRange'), opt.FrameRange = [ceil(opt.TOIsize/2), ceil(opt.TOIsize/2)+1];end
    if ~isfield(opt, 'SaveData'), opt.SaveData = 1; end
    if ~isfield(opt, 'ProgressBar'), opt.ProgressBar = 1; end
    if ~isfield(opt, 'Parallel'), opt.Parallel = 1; end
    if ~isfield(opt, 'minR2'), opt.minR2 = 0.99; end
    if ~isfield(opt, 'InferNaNs'), opt.InferNaNs =1; end
    if ~isfield(opt, 'Interpolate'), opt.Interpolate =1; end

    % shift of ROI's amount & define the position of ROI-TOI
    opt.fracROIshift = opt.ROIsize/opt.ROIshift;
    opt.fracTOIshift = opt.TOIsize/opt.TOIshift;

    along_y = floor((size(stack1,2)-opt.ROIsize)/opt.ROIshift)+1;
    along_x = floor((size(stack1,1)-opt.ROIsize)/opt.ROIshift)+1;
    along_t = floor((size(stack1,3)-opt.TOIsize)/opt.TOIshift)+1;

    if ~opt.AllFrames
        if isempty(opt.FrameRange)
            opt.AllFrames = 1;
            kRange = 1:along_t;
            opt.FrameRange = kRange + floor(opt.TOIsize/2);
        else
            if mod(opt.TOIsize, 2) ==  0 % TOI size is even
                minK = 1;
            else
                minK = 0;
            end
            kRange = opt.FrameRange - floor(opt.TOIsize/2) + minK; % Vector field frame
            if length(kRange) == 1; opt.Parallel = 0; end % No use running in parallel if its only one frame
        end
    else
        kRange = 1:along_t;
        opt.FrameRange = kRange + floor(opt.TOIsize/2);        
    end
    opt.kRange = kRange;

    opt.position_x  =   zeros(1,along_x);
    opt.position_y  =   zeros(1,along_y);
    opt.position_t  =   zeros(1,along_t);

    % defining all of the ROI centers positons in x and y
    for i=1:along_x
        opt.position_x(i) = 1/2 + (i-1)/opt.fracROIshift*opt.ROIsize + opt.ROIsize/2;
    end
    opt.position_x = repmat(opt.position_x,along_y,1);
    opt.position_x = opt.position_x';

    for j=1:along_y
        opt.position_y(j) = 1/2 + (j-1)/opt.fracROIshift*opt.ROIsize + opt.ROIsize/2;
    end
    opt.position_y = repmat(opt.position_y,along_x,1);


%     if size(opt.maskCell,3) == 1  % Static mask
        % %pre-process to determine if you want to analyse the whole FOV or select a
        % %polygon to analyze...in case of polygon, it will speed up the whole
        % %calculations because less ROI-TOI to consider...
        % imagesc(mean(series,3)) 
        % prompt = {'want to select the polygon or use whole FOV?'};
        %     dlg_title = 'y or n';
        %     num_lines = 1;
        %     def = {'y'};
        %     answer = inputdlg(prompt,dlg_title,num_lines,def);
        % if strcmp(answer,'y')
        %     maskCell=[];
        %     while isempty(maskCell)
        %         [x,y,maskCell,poly1.x,poly1.y]=roipoly;
        %     end
        % elseif strcmp(answer,'n')
        %     maskCell=ones(size(series,1),size(series,2));
        % end
        % close
        % opt.maskCell = maskCell;

        % %% Old way: if at least one pixel of the ROI is outside of the mask, then
        % % don't calculate that vector.
        % vectorPositions = zeros(along_x,along_y);
        % for w=1:along_x
        %     for s=1:along_y
        %         indx = 1+(w-1)/opt.fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((w-1)/opt.fracROIshift)));
        %         indy = 1+(s-1)/opt.fracROIshift*opt.ROIsize:(opt.ROIsize*(1+((s-1)/opt.fracROIshift))); 
        %         indx = floor(indx);
        %         indy = floor(indy);
        %         if(opt.maskCell(indx,indy)==1)
        %             vectorPositions(w,s)=1;
        %         end
        %         if(opt.dynamicMask(indx,indy)==1)
        %             maskedPositions(w,s)=1;
        %         end        
        %     end
        % end
        % opt.vectorPositions = vectorPositions;

        % New way: plot all those vectors whose positions fall within the mask,
        % even if the ROI extends further beyond the mask
        vectorPositions = zeros(along_x, along_y);
        maskedPositions = zeros(along_x, along_y, along_t);
        for w = 1:along_x
            for s = 1:along_y
                indx = opt.position_x(w,1);
                indy = opt.position_y(1,s); 
                indx = floor(indx);
                indy = floor(indy);
                if(opt.maskCell(indx,indy)==1)
                    vectorPositions(w,s)=1;
                end
                
                for t = 1: along_t                    
                    if(opt.dynamicMask(indx, indy, t)==1)
                        maskedPositions(w,s, t)=1;
                    end
                end
            end
        end
        opt.vectorPositions = vectorPositions;
        opt.maskedPositions = logical(maskedPositions);

    % %indexes of ROI within the mask of interest (polygon)
    [iroi, jroi] = find(vectorPositions);
    velocityMap = cell(1, along_t);
        

    %% Immobile filtering
    if strcmp(opt.filtering,'FourierWhole')
        stack1 = immfilter(stack1,'F',1); 
        stack2 = immfilter(stack2,'F',1); 
    elseif strcmp(opt.filtering,'MovingAverage')
        stack1 = immfilter(stack1,opt.MoveAverage); 
        stack2 = immfilter(stack2,opt.MoveAverage); 
    elseif strcmp(opt.filtering,'butterIIR')
        % set the cutoff frequency to be the 2*sampling frequency / size of the series 
        % this way we will certainly remove the immobile component
        stack1 = butterIIR(stack1, 2/(opt.timeFrame*size(stack1,3)),'none');
        stack2 = butterIIR(stack2, 2/(opt.timeFrame*size(stack1,3)),'none');
    elseif strcmp(opt.filtering,'none')
    else
        error('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
    end
    
    %% Progress bar
    if opt.ProgressBar
        if opt.Parallel
            try
                h = ICS_ExplorerUIFigure;
            catch
                h = uifigure;                
            end
            d = uiprogressdlg(h,"Title", "Calculating STICS in Parallel...", "Message",...
                {['Calculating ', num2str(length(kRange)),' vector fields in parallel with ', num2str(length(iroi)), ' vectors per field'],...
                ['for a total of ', num2str(length(iroi)*length(kRange)), ' vectors. Please wait.']});
            d.Indeterminate = 1; 
            drawnow;
        else        
            h = uifigure;
            d = uiprogressdlg(h, 'Title', 'Calculating STICS...','Message',...
                                {['Vector Field #0, (0/', num2str(length(kRange)),...
                                ' = ', num2str(round(100*(0)/length(kRange),1)),'%).'],...
                                ['Vector #0, (0/',num2str(length(iroi)*length(kRange)), ' = ',...
                                num2str(round(100*(((0)*length(iroi)))/(length(iroi)*length(kRange)),1)),'%)']},...
                                'Indeterminate', 0);
                            d.ShowPercentage = 1;
        end
    end
    tEnd = 0;
    
    %% Define a nested function to generate each vector field
    function VFstruct = generateVectorField(k, cfn)

        TOIFOV1 = stack1(:,:,floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift))));
        TOIFOV2 = stack2(:,:,floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift))));
        % define the average FOV TOI image that will be used in the display of the vector maps
%         VFstruct.data_TOImean = mean(TOIFOV1, 3); % -> Removed so the
%         output is less heavy and the calculation requires less memory.

        VFstruct.vx = zeros(size(opt.position_x));
        VFstruct.vy = zeros(size(opt.position_x));

        VFstruct.R2x = zeros(size(opt.position_x));
        VFstruct.R2y = zeros(size(opt.position_y));        

        StartTime = now; tStart = tic; % 
        fprintf('Analyzing TOI %i of %i ____ date and time is %s',k,length(kRange),datestr(now)); 

        VFstruct.sigFitRatio = zeros(along_x, along_y);
    %     VFstruct.coeffGtime = cell(along_x, along_y);
    %     VFstruct.Vx = zeros(along_x, along_y, 2);
    %     VFstruct.Vy = zeros(along_x, along_y, 2);

        for u=1:length(iroi) 
            
            % Update progress bars
            if opt.ProgressBar
                if ishandle(h)
                    if ~opt.Parallel
                        
                        d.Value = (((find(kRange == k)-1)*length(iroi)) + u)/(length(iroi)*length(kRange));
                        d.Message = {['Working on vector field #', num2str(k),', (' num2str(find(kRange == k)-1), '/', num2str(length(kRange)),...
                                ' = ', num2str(round(100*(find(kRange == k)-1)/length(kRange),1)),'% done).'],...
                                ['Vector # ', num2str(u), ' / ', num2str(length(iroi)), ' (',...
                                num2str(round(100*(u)/(length(iroi)),1)),'% done).'],...
                                ['Last vector field took ', num2str(floor(tEnd/60)), ' minutes and ', num2str(floor(rem(tEnd,60))), ' seconds.'],
                                ['Start time: ' datestr(StartTime)]};    
%                                  For overall progress:
%                                 ['Vector # ', num2str((((find(kRange == k))-1)*length(iroi)) + u), ' (',...
%                                 num2str((((find(kRange == k))-1)*length(iroi)) + u),'/',...
%                                 num2str(length(iroi)*length(kRange)), ' = ',...
%                                 num2str(round(100*(((find(kRange == k)-1)*length(iroi)) + u)/(length(iroi)*length(kRange)),1)),'% done).'],...                                
%                             ['At that rate, we estimate that the analysis would be done by ' datestr(StartTime + (length(kRange)*tEnd))]
                        drawnow;

                    end
                end
            end        

            % Accessing the data of only the ROIs that are defined within the selected polygon
            i=iroi(u);
            j=jroi(u);           
            
            % Define the data region (x,y,t) on which STICS will be applied
            FOV1 = TOIFOV1(floor(1+(i-1)/opt.fracROIshift*opt.ROIsize):floor(opt.ROIsize*(1+((i-1)/opt.fracROIshift))), floor(1+(j-1)/opt.fracROIshift*opt.ROIsize): floor(opt.ROIsize*(1+((j-1)/opt.fracROIshift))),:);
            FOV2 = TOIFOV2(floor(1+(i-1)/opt.fracROIshift*opt.ROIsize):floor(opt.ROIsize*(1+((i-1)/opt.fracROIshift))), floor(1+(j-1)/opt.fracROIshift*opt.ROIsize): floor(opt.ROIsize*(1+((j-1)/opt.fracROIshift))),:);

            % Apply regular cross-corr STICS
            switch opt.STICCS
                case 1
                    [corrfn12, corrfn21, corrfn1, corrfn2] = stics(FOV1, FOV2, opt.tauLimit);
                    corrfns = {corrfn1, corrfn2, corrfn12, corrfn21};
                    corrfn = corrfns{1,cfn};
                case 0
                    [corrfn1] = stics(FOV1, opt.tauLimit);
                    corrfn = corrfn1; clear corrfn1;
            end

            % Check for significance of global maximum
            upperTauLimit = min(opt.tauLimit, size(FOV1,3)); % Since size(stack1)==size(stack1), size(FOV1)==size(FOV2)
            if opt.CheckSignificance
                for  tau = 0:upperTauLimit-1                                      
                    if ~correlationSignificance(corrfn(:,:,tau+1))
                        corrfn = corrfn(:,:,1:(end-1)); % cut off the "bad" lag(s)
%                         corrfns{1,cfn} = corrfn;                                                    
                        break
                    end
                end            
            end            
            

            while VFstruct.R2x(i,j) < opt.minR2 || VFstruct.R2y(i,j) < opt.minR2

                % Fit vxtau and vytau vs tau linearly to extract vx and vy
                if size(corrfn,3)>1 && size(corrfn,1)==opt.ROIsize && size(corrfn,2)==opt.ROIsize % second two arguments added MM because of error when corrfn size is strange
                    % Apply symmetric gaussian fit to data
                    [coeffGtime] = gaussfit(corrfn,'time',opt.pixelSize,'n',opt.fitRadius);
    %                 [coeffGrotated] = gaussfit(corrfn(:,:,1),'rotated',opt.pixelSize,'n',opt.fitRadius);


    %% Quality Control
                    % Discard correlation function fits if they have:
                    % negative amplitudes or...
                    % large X width or...
                    % large Y width...
                    % or an amplitude decay beyond 1/4.
                    % or an amplitude greater than the first.
                    
                    % Remove these fits and their corresponding correlation functions:
                    [coeffGtime, corrfn] = qualityControlCoeffGfit(coeffGtime, corrfn, opt);


    %% Linear fit model
                    % If there are more than 2 fits left:
                    if size(corrfn,3) >= 2 && size(coeffGtime,1) >= 2
                        % Calculate linear model fit
                        [VFstruct.Px(u,2:-1:1), VFstruct.Sx(u)] = polyfit(opt.timeFrame*(1:size(corrfn,3)), squeeze(coeffGtime(1:size(corrfn,3),5))', 1);
                        [VFstruct.Py(u,2:-1:1), VFstruct.Sy(u)] = polyfit(opt.timeFrame*(1:size(corrfn,3)), squeeze(coeffGtime(1:size(corrfn,3),6))', 1);
                    else
                        VFstruct.Px(u,2:-1:1) = NaN; 
                        VFstruct.Py(u,2:-1:1) = NaN; 
                    end
                else
                    VFstruct.Px(u,2:-1:1) = NaN; 
                    VFstruct.Py(u,2:-1:1) = NaN; 
                end

                if ~isnan(VFstruct.Px(u,2)) && ~isnan(VFstruct.Py(u,2))
                    if size(corrfn,3) >= 3 && size(coeffGtime,1) >= 3                
                        [Ex, ~] = polyval(VFstruct.Px(u,2:-1:1), opt.timeFrame*(1:size(corrfn,3)), VFstruct.Sx(u));
                        xresid = squeeze(coeffGtime(1:size(corrfn,3),5))' - Ex;
                        SSxresid = sum(xresid.^2);
                        SSxtotal = (length(squeeze(coeffGtime(1:size(corrfn,3),5))')-1 * var(squeeze(coeffGtime(1:size(corrfn,3),5))'));
                        xrsq = 1 - SSxresid/SSxtotal;
                        VFstruct.R2x(i,j) = xrsq;
    %                     VFstruct.Sx.DeltaX(u) = deltax;
                        VFstruct.SSrx(i,j) = SSxresid;

                        [Ey, ~] = polyval(VFstruct.Py(u,2:-1:1), opt.timeFrame*(1:size(corrfn,3)), VFstruct.Sy(u));
                        yresid = squeeze(coeffGtime(1:size(corrfn,3),6))' - Ey;
                        SSyresid = sum(yresid.^2);
                        SSytotal = (length(squeeze(coeffGtime(1:size(corrfn,3),6))')-1 * var(squeeze(coeffGtime(1:size(corrfn,3),6))'));
                        yrsq = 1 - SSyresid/SSytotal;
                        VFstruct.R2y(i,j) = yrsq;
    %                     VFstruct.Sy.DeltaY(u) = deltay;
                        VFstruct.SSry(i,j) = SSyresid;
                    elseif size(corrfn,3) == 2 && size(coeffGtime,1) == 2
                        VFstruct.R2x(i,j) = 1;
                        VFstruct.SSrx(i,j) = 0;
                        VFstruct.R2y(i,j) = 1;                        
                        VFstruct.SSry(i,j) = 0;                        
                    end
                else
                    VFstruct.R2x(i,j) = NaN;
    %                 VFstruct.Sx.DeltaX(u) = NaN;
                    VFstruct.SSrx(i,j) = NaN;

                    VFstruct.R2y(i,j) = NaN;
    %                 VFstruct.Sy.DeltaY(u) = NaN;
                    VFstruct.SSry(i,j) = NaN;                                             
                end

                if VFstruct.R2x(i,j) < opt.minR2 || VFstruct.R2y(i,j) < opt.minR2
                    corrfn = corrfn(:,:,1:end-1);              
                end
            end
            
            VFstruct.vx(i,j) = VFstruct.Px(u,2);
            VFstruct.vy(i,j) = VFstruct.Py(u,2);              
            
            VFstruct.sigFitRatio(i,j) = size(corrfn,3)/size(FOV1,3);
            
        end % end for indx

        VFstruct.goodVectors = logical(~isnan(VFstruct.vx) & ~isnan(VFstruct.vy)); 
        
        if opt.InferNaNs
            VFstruct.WasNaN = (isnan(VFstruct.vx) | isnan(VFstruct.vy));            
            VFstruct.vx = inpaint_nans(VFstruct.vx);
            VFstruct.vx(VFstruct.vx == 0) = NaN;
            VFstruct.vy = inpaint_nans(VFstruct.vy);     
            VFstruct.vy(VFstruct.vy == 0) = NaN;
        end                
        
        VFstruct.maskedVectors = logical(VFstruct.goodVectors.*opt.maskedPositions(:,:,k));
        
        if ~strcmp(opt.VelUnits, 'um/min') || ~isfield(VFstruct, "VelMap")
            switch opt.TimeUnits
                case 'sec'
                    VFstruct.VelMap = sqrt(VFstruct.vx.^2 + VFstruct.vy.^2).*60;
                    VFstruct.vx = VFstruct.vx.*60;
                    VFstruct.vy = VFstruct.vy.*60;                    
                case 'min'
                    VFstruct.VelMap = sqrt(VFstruct.vx.^2 + VFstruct.vy.^2);
            end
        end                     
                                 
        
        if opt.Interpolate
            VFstruct.iVx = zeros(size(opt.maskCell));
            VxMap = imresize(VFstruct.vx, opt.ROIshift, 'bilinear');
            VFstruct.iVx(ceil(opt.ROIsize/2):floor(opt.ROIsize/2)+size(VxMap,1)-1,...
                ceil(opt.ROIsize/2):floor(opt.ROIsize/2)+size(VxMap,2)-1) = VxMap;
            
            VFstruct.iVy = zeros(size(opt.maskCell));            
            VyMap = imresize(VFstruct.vy, opt.ROIshift, 'bilinear');
            VFstruct.iVy(ceil(opt.ROIsize/2):floor(opt.ROIsize/2)+size(VyMap,1)-1,...
                ceil(opt.ROIsize/2):floor(opt.ROIsize/2)+size(VyMap,2)-1) = VyMap;
            
            VFstruct.iVelMap = sqrt(VFstruct.iVx.^2 + VFstruct.iVy.^2);
            
        end

        tEnd = toc(tStart);
        fprintf(' ____ TOI finished in %d minutes and %f seconds\n',floor(tEnd/60),floor(rem(tEnd,60))); % added MM

    end
    
    % Create a handle for the function to use it in parallel.
    GVF = @generateVectorField;
    
    %% Run FOR or PARFOR loop across the time interval
    % Define the number of outputs:
    switch opt.STICCS
        case 1
            nOutputs = 4;
        case 0
            nOutputs = 1;
    end
    
    % Define the position of TOI centers in time:
    for k = kRange
        opt.position_t(k) = 1/2 + (k-1)/opt.fracTOIshift*opt.TOIsize + opt.TOIsize/2;     
    end      
    
    if opt.Parallel                      
        for cfn = 1:nOutputs
            parfor k = kRange                                   
                velocityMap{cfn,k} = feval(GVF, k, cfn); %#ok<FVAL>
            end
        end
    else    
        for cfn = 1:nOutputs
            for k = kRange             
                velocityMap{cfn,k} = generateVectorField(k, cfn);
            end
        end
    end
    opt.VelUnits = 'um/min';    
    
    %% Save data
    if opt.ProgressBar
        d.Indeterminate = 1;
        d.Message = 'Saving data...';
        drawnow;
    end  

    if opt.SaveData
        % output data to file
        save([opt.path opt.outputName 'VM.mat'],'velocityMap',  '-v7.3');
        save([opt.path opt.outputName 'VM.mat'],'opt','-append');
    end

    if opt.ProgressBar
        if ishandle(h)
            delete(h)
        end
    end
end 

