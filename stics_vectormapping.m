function [velocityMap, position_x, position_y, position_t, opt] = stics_vectormapping(stack, opt)
%this functions takes the timeseries and calculates for defined regions of
%size 'ROIsize' and 'TOIsize' stics correlation function and fits symmetric
%Gaussian to extract the flow in x and y directions...
%Developed by Elvis Pandzic 2012-2014 under supervision of 
%Prof. Paul W. Wiseman @ McGill University. 
% correspondence addressed to : e.pandzic@unsw.edu.au


    % Look at opt structure and either assign default value, or use passed
    % value
    if nargin==1
        opt.pixelSize = 0.1; % just initialize it to something so the code below doesn't throw an error
    end
    if ~isfield(opt, 'pixelSize'), opt.pixelSize = 1; end %pixel size in um
    if ~isfield(opt, 'timeFrame'), opt.timeFrame = 1; end %time delay (s) between subsequent frames
    if ~isfield(opt, 'TimeUnits'), opt.TimeUnits = 'sec'; end
    if ~isfield(opt, 'tauLimit'), opt.tauLimit = size(stack,3); end % what is highest time lag to compute stics of TO
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
    if ~isfield(opt, 'maskCell'), opt.maskCell = ones(size(stack, 1), size(stack, 2)); end

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

    % shift of ROI's amount & define the position of ROI-TOI
    opt.fracROIshift = opt.ROIsize/opt.ROIshift;
    opt.fracTOIshift = opt.TOIsize/opt.TOIshift;

    along_y = floor((size(stack,2)-opt.ROIsize)/opt.ROIshift)+1;
    along_x = floor((size(stack,1)-opt.ROIsize)/opt.ROIshift)+1;
    along_t = floor((size(stack,3)-opt.TOIsize)/opt.TOIshift)+1;

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
        end
    else
        kRange = 1:along_t;
        opt.FrameRange = kRange + floor(opt.TOIsize/2);        
    end

    position_x  =   zeros(1,along_x);
    position_y  =   zeros(1,along_y);
    position_t  =   zeros(1,along_t);

    % defining all of the ROI centers positons in x and y
    for i=1:along_x
        position_x(i) = 1/2 + (i-1)/opt.fracROIshift*opt.ROIsize + opt.ROIsize/2;
    end
    position_x=repmat(position_x,along_y,1);
    position_x=position_x';

    for j=1:along_y
        position_y(j) = 1/2 + (j-1)/opt.fracROIshift*opt.ROIsize + opt.ROIsize/2;
    end
    position_y=repmat(position_y,along_x,1);


    if size(opt.maskCell,3) == 1  % Static mask
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
        %     end
        % end
        % opt.vectorPositions = vectorPositions;

        % New way: plot all those vectors whose positions fall within the mask,
        % even if the ROI extends further beyond the mask
        vectorPositions = zeros(along_x,along_y);
        for w=1:along_x
            for s=1:along_y
                indx = position_x(w,1);
                indy = position_y(1,s); 
                indx = floor(indx);
                indy = floor(indy);
                if(opt.maskCell(indx,indy)==1)
                    vectorPositions(w,s)=1;
                end
            end
        end
        opt.vectorPositions = vectorPositions;
    else
        dynamicMask = 1;
    end


    %immobile filter data if necessary
    if strcmp(opt.filtering,'FourierWhole')
        stack = immfilter(stack,'F',1); 
    elseif strcmp(opt.filtering,'MovingAverage')
        stack = immfilter(stack,opt.MoveAverage); 
    elseif strcmp(opt.filtering,'butterIIR')
        %set the cutoff frequency to be the 2*sampling frequency / size of the series 
        % this way we will certainly remove the immobile component
        stack = butterIIR(stack, 2/(opt.timeFrame*size(stack,3)),'none');
    elseif strcmp(opt.filtering,'none')
    else
        error('Please select a valid filter (''FourierWhole'' ,''MovingAverage'' , ''butterIIR'' or ''none'')');
    end

    % %indexes of ROI within the mask of interest (polygon)
    [iroi, jroi] = find(vectorPositions);
    velocityMap = cell(1, along_t);

    if opt.ProgressBar
        if opt.Parallel
            h = uifigure;
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

    
    %% Define a nested function to generate each vector field
    function VFstruct = generateVectorField(k)

        TOIFOV = stack(:,:,floor(1+(k-1)/opt.fracTOIshift*opt.TOIsize):floor(opt.TOIsize*(1+((k-1)/opt.fracTOIshift))));
        % define the average FOV TOI image that will be used in the display of the vector maps
        VFstruct.data_TOImean = mean(TOIFOV, 3);



        tStart = tic; % 
        fprintf('analyzing TOI %i of %i ____ date and time is %s',k,length(kRange),datestr(now)); 

        VFstruct.sigFitRatio = zeros(along_x, along_y);
    %     VFstruct.coeffGtime = cell(along_x, along_y);
    %     VFstruct.Vx = zeros(along_x, along_y, 2);
    %     VFstruct.Vy = zeros(along_x, along_y, 2);
    %     VFstruct.VelMap = zeros(along_x, along_y);

        for u=1:length(iroi) 
            
            % Update progress bars
            if opt.ProgressBar
                if ishandle(h)
                    if ~opt.Parallel
                        waitbar((((find(kRange == k)-1)*length(iroi)) + u)/(length(iroi)*length(kRange)), h,...
                            {['Vector Field #', num2str(k),', (' num2str(find(kRange == k)-1), '/', num2str(length(kRange)),...
                                ' = ', num2str(round(100*(find(kRange == k)-1)/length(kRange),1)),'%).'],
                                ['Vector #', num2str((((find(kRange == k))-1)*length(iroi)) + u), ' (',...
                                num2str((((find(kRange == k))-1)*length(iroi)) + u),'/',...
                                num2str(length(iroi)*length(kRange)), ' = ',...
                                num2str(round(100*(((find(kRange == k)-1)*length(iroi)) + u)/(length(iroi)*length(kRange)),1)),'%)']});
                    end
                end
            end        

            % Accessing the data of only the ROIs that are defined within the selected polygon
            i=iroi(u);
            j=jroi(u);
            
            % Define the data region (x,y,t) on which STICS will be applied
            regionanalyse = TOIFOV(floor(1+(i-1)/opt.fracROIshift*opt.ROIsize):floor(opt.ROIsize*(1+((i-1)/opt.fracROIshift))), floor(1+(j-1)/opt.fracROIshift*opt.ROIsize): floor(opt.ROIsize*(1+((j-1)/opt.fracROIshift))),:);

            % Apply regular cross-corr STICS
            [corrfn] = stics(regionanalyse, opt.tauLimit);

            % Check for significance of global maximum
            upperTauLimit = min(opt.tauLimit, size(regionanalyse,3));
            if opt.CheckSignificance
                for  tau = 0:upperTauLimit-1
                    if ~correlationSignificance(corrfn(:,:,tau+1))
                        corrfn = corrfn(:,:,1:(end-1)); % cut off the "bad" lag(s)
                        break
                    end            
                end            
            end            

            % Fit vxtau and vytau vs tau linearly to extract vx and vy
            if size(corrfn,3)>1 && size(corrfn,1)==opt.ROIsize && size(corrfn,2)==opt.ROIsize % second two arguments added MM because of error when corrfn size is strange
                % Apply symmetric gaussian fit to data
                [coeffGtime] = gaussfit(corrfn,'time',opt.pixelSize,'n',opt.fitRadius);
%                 [coeffGrotated] = gaussfit(corrfn(:,:,1),'rotated',opt.pixelSize,'n',opt.fitRadius);
%                 VFstruct.coeffGtime{i,j} = coeffGtime;
%                 VFstruct.coeffGrotated{i,j} = coeffGrotated;


%% Quality Control
                % Discard correlation function fits if they have:
                cutOff = min([find(coeffGtime(:,1)<0, 1, 'first'),... % negative amplitudes or...
                    find(coeffGtime(:,2)>opt.maxHalfWidth,1,'first'),... % large X width or...
                    find(coeffGtime(:,3)>opt.maxHalfWidth,1,'first'),... % large Y width...
                    find(coeffGtime(:,1)<coeffGtime(1,1)/4, 1, 'first'),... % or an amplitude decay beyond 1/4.
                    find(coeffGtime(:,1)>coeffGtime(1,1).*1.1, 1, 'first')]); % or an amplitude greater than the first.

                % Remove these fits and their corresponding correlation
                % functions:
                if cutOff >= 1
                    coeffGtime = coeffGtime(1:cutOff-1,:);
                    corrfn = corrfn(:,:,1:cutOff-1);
                end
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
            VFstruct.sigFitRatio(i,j) = size(corrfn,3)/size(regionanalyse,3);
        end % end for indx

        %% Reorganize Px and Py into vx and vy:
        vx = nan(size(position_x));
        vy = nan(size(position_y));
        
        for u=1:length(iroi)
            i=iroi(u);
            j=jroi(u);
            vy(i,j)=squeeze(VFstruct.Py(u,2));
            vx(i,j)=squeeze(VFstruct.Px(u,2));
        end
        
        VFstruct.vx=vx;
        VFstruct.vy=vy;

        VFstruct.goodVectors = (~isnan(vx) & ~isnan(vy)); 

        tEnd = toc(tStart);
        fprintf(' ____ TOI finished in %d minutes and %f seconds\n',floor(tEnd/60),round(rem(tEnd,60))); % added MM

    end
    
    % Create a handle for the function to use it in parallel.
    GVF = @generateVectorField;
    
    %% Run FOR or PARFOR loop across the time interval
    if opt.Parallel
        
        for k = kRange%1:along_t % loop along time                        
            % position of TOI centers in time
            position_t(k) = 1/2 + (k-1)/opt.fracTOIshift*opt.TOIsize + opt.TOIsize/2;
            % define whole FOV TOI            
        end                
        
        parfor k = kRange                                   
            velocityMap{k} = feval(GVF, k); %#ok<FVAL>
        end
        
    else    
        for k = kRange%1:along_t % loop along time                        
            % position of TOI centers in time
            position_t(k) = 1/2 + (k-1)/opt.fracTOIshift*opt.TOIsize + opt.TOIsize/2;
            % define whole FOV TOI                        
            velocityMap{k} = generateVectorField(k);
        end
    end
    
    %% Save data
    if opt.ProgressBar
        d.Indeterminate = 1;
        d.Message = 'Saving data...';
        drawnow;
    end  

    if opt.SaveData
        % output data to file
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'velocityMap',  '-v7.3');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_x','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_y','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_t','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'opt','-append');
    end

    if opt.ProgressBar
        if ishandle(h)
            delete(h)
        end
    end
end 

