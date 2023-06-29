function [velocityMap, opt] = convertVelocitiesFromPxPy(velocityMap, position_x, position_y, position_t, opt)

    if ~isfield(opt, 'TimeUnits'), opt.TimeUnits = 'sec'; end
    if ~isfield(opt, 'VelUnits'), opt.VelUnits = 'um/sec'; end
    if ~isfield(opt, 'AllFrames'), opt.AllFrames = 1; end
    if ~isfield(opt, 'FrameRange'), opt.FrameRange = [ceil(opt.TOIsize/2), size(velocityMap, 2) + floor(opt.TOIsize/2)];end%[ceil(opt.TOIsize/2), ceil(opt.TOIsize/2)+1];end    
    if ~strcmp(opt.VelUnits, 'um/min') || ~isfield(velocityMap{1}, "VelMap")
        if ~opt.AllFrames
            if isempty(opt.FrameRange)
                opt.AllFrames = 1;
            else
                if mod(opt.TOIsize, 2) ==  0 % TOI size is even
                    minK = 1;
                else
                    minK = 0;
                end
                kRange = opt.FrameRange - floor(opt.TOIsize/2) + minK; % Vector field frame
            end
        else
            kRange = 1:length(velocityMap);
        end

        [iroi, jroi] = find(opt.vectorPositions);    
        vx=nan(size(position_x));
        vy=nan(size(position_y));

        for k=kRange % loop along time
            if ~isfield(velocityMap{k}, 'vx') || ~isfield(velocityMap{k}, 'vy')     
                for u=1:length(iroi)
                    i=iroi(u);
                    j=jroi(u);
                    if isfield(velocityMap{k}, "Px")
                        vy(i,j)=squeeze(velocityMap{k}.Py(u,2));
                        vx(i,j)=squeeze(velocityMap{k}.Px(u,2));
                    end
                end
                velocityMap{k}.vx=vx;
                velocityMap{k}.vy=vy;        
            end

            switch opt.TimeUnits
                case 'sec'
                    velocityMap{k}.VelMap = sqrt(velocityMap{k}.vx.^2 + velocityMap{k}.vy.^2).*60;
                    velocityMap{k}.vx = velocityMap{k}.vx.*60;
                    velocityMap{k}.vy = velocityMap{k}.vy.*60;
                case 'min'
                    velocityMap{k}.VelMap = sqrt(velocityMap{k}.vx.^2 + velocityMap{k}.vy.^2);
            end
        end
    end

    opt.VelUnits = 'um/min';
    if ~isfield(opt, 'SaveData'), opt.SaveData = 0; end

    if opt.SaveData
        % output data to file
        save([opt.path opt.outputName 'VM.mat'],'velocityMap',  '-v7.3');
        save([opt.path opt.outputName 'VM.mat'],'position_x','-append');
        save([opt.path opt.outputName 'VM.mat'],'position_y','-append');
        save([opt.path opt.outputName 'VM.mat'],'position_t','-append');
        save([opt.path opt.outputName 'VM.mat'],'opt','-append');
    end        

end