function [velocityMap, opt] = setupForPlotting(velocityMap,position_t,position_x,position_y,opt)



    for i = 1:length(velocityMap)

        velocityMap{i}.goodVectors(end) = 1; % Force last vector to be valid for auto scale

        vy = squeeze(velocityMap{i}.vy);
        vyall(i,:) = vy(1:size(vy,1)*size(vy,2));
        vx = squeeze(velocityMap{i}.vx);
        vxall(i,:) = vx(1:size(vx,1)*size(vx,2));

        velocityMap{i}.goodVectors = reshape(velocityMap{i}.goodVectors,...
            size(squeeze(position_x,1)), size(squeeze(position_y,2)));
        tempMap = velocityMap{i}.VelMap.*velocityMap{i}.goodVectors;
        AllVelMap(:,i) = tempMap(1:size(tempMap,1)*size(tempMap,2));
    end

    switch opt.VelScaleRange
        case 'Global'
            upperVLimit = round(prctile(AllVelMap, 99.75, "all"),2);

            % Sneak the all times maximum into each vector field so that the auto scale
            % quiver is adjusted appropriately to compare colors accross frames.
            vxall(:,end) = -sqrt(upperVLimit);
            vyall(:,end) = -sqrt(upperVLimit);
        case 'Local'

            upperVLimit = round(prctile(velocityMap{k}.VelMap.*velocityMap{k}.goodVectors, 99.75, "all"),2);
            % OR
            % upperVLimit = round(max(velocityMap{1}.VelMap,[],"all"),2);
    end

    upperVLimit = min([upperVLimit, opt.maxV]);
    
    indgoodvect = find(velocityMap{1}.goodVectors);
    overallgoodvectors = indgoodvect;
    allgood = [];

    for i = 1:length(velocityMap)
        % all vector that magnitudes per minutes smaller than certain thesholdV
        badVectors = find(velocityMap{i}.VelMap>upperVLimit);

        velocityMap{i}.plottedVectors = setdiff(indgoodvect,badVectors);
        allgood = union(allgood,velocityMap{i}.plottedVectors);
        overallgoodvectors = intersect(overallgoodvectors,velocityMap{i}.plottedVectors);

        velocityMap{i}.minVel = min(velocityMap{i}.VelMap(velocityMap{i}.plottedVectors), [],"omitnan","all");
        velocityMap{i}.maxVel = max(velocityMap{i}.VelMap(velocityMap{i}.plottedVectors), [],"omitnan","all");

    end

    minGlobalVel = min(min(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");
    maxGlobalVel = max(max(AllVelMap(overallgoodvectors,:), [], "omitnan"), [], "omitnan");

    opt.GlobalVelRange = [minGlobalVel, upperVLimit, maxGlobalVel];
    
    if opt.SaveData
        % output data to file
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'velocityMap',  '-v7.3');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_x','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_y','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'position_t','-append');
        save([opt.path 'VelocityMap' opt.outputName '.mat'],'opt','-append');
    end    
end