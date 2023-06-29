%% Generate significance mean and min maps from several velocity maps files
% Rodrigo Migueles April 2nd 2021
regions = 1:6; ROIsizes = [8, 16, 32, 64]; TOIsizes = [3, 5, 10];
sigContainer = cell(length(ROIsizes)*length(TOIsizes), length(regions)+2);
nfileIx = 0;
for region = 1:6
    disp(['Region #', num2str(region)]);
    [fileList, filePath] = uigetfile(".mat","MultiSelect","on");
    if ischar(fileList) % Only one file selected
        fileList = {fileList};
    end

    % figure;
    for fileIx = 1:length(fileList)
        load([filePath, fileList{fileIx}]); 
        fileName = fileList{fileIx}; fileName = fileName(1:end-4);

        sigContainer{fileIx, 1} = opt.ROIsize; sigContainer{fileIx, 2} = opt.TOIsize;
        sigContainer{fileIx, region+2} = zeros(size(position_x,1), size(position_y,2), length(velocityMap));
        for f = 1:length(velocityMap)
            sigContainer{fileIx, region+2}(:,:,f) = velocityMap{f}.sigFitRatio;
        end
    %     sigMeanMap = mean(sigContainer{fileIx}, 3);
    %     sigMinMap = min(sigContainer{fileIx}, [], 3);

    %     imagesc(sigContainer(:,:,2), [0,1]);axis image; colorbar; ax = gca;
    %     exportgraphics(ax, [filePath, [fileName, '_sigK2Map.png']]);

    %     subplot(5, 4, fileIx)
    %     title(['TOI size = ', num2str(opt.TOIsize)]);
    %     imagesc(sigMinMap, [0, 1]); axis image; 
    %     ylabel(['ROI size = ', num2str(opt.ROIsize)]);
    %     xlabel(['TOI size = ', num2str(opt.TOIsize)]);
    %     ax = gca; exportgraphics(ax, [filePath, [fileName, '_sigMeanMap.png']]);

    %     imagesc(sigMinMap, [0, 1]); axis image; colorbar; ax = gca;
    %     exportgraphics(ax, [filePath, [fileName, '_sigMinMap.png']]);
    end
    nfileIx = nfileIx + length(fileList);
end

%%
% 

for r = regions
    i = 1;
    for ROIsize = 1:4
        for TOIsize = 1:3
            MeansOfMeans(ROIsize, TOIsize, r) = mean(sigContainer{i, r+2}, 'all');
            MeansOfMins(ROIsize, TOIsize, r) =  mean(min(sigContainer{i, r+2}, [], 3), 'all');
            MeansOfVars(ROIsize, TOIsize, r) =  mean(var(sigContainer{i, r+2}, 0, 3), 'all');
            i = i + 1;
        end
    end
end
%%
for p = 1:3
    switch p
        case 1
            mapToUse = MeansOfMeans;
        case 2
            mapToUse = MeansOfMins;
        case 3
            mapToUse = MeansOfVars;
    end
    
    figure;
    for r = regions
        subplot(2,3,r)
        imagesc(mapToUse(:,:,r)); axis square; colorbar
        title(['Region ', num2str(r)]);
        xlabel('TOI size'); xticklabels({'3', '5', '10'});
        ylabel('ROI size'); yticklabels({'8', '16', '32', '64'});
    end
    switch p
        case 1
            sgtitle('Mean of mean significance ratio maps');
        case 2
            sgtitle('Mean of minimum significance ratio maps');
        case 3
            sgtitle('Mean of variance significance ratio maps');
    end
%     a = gcf;
%     exportgraphics(a,
end

% sgtitle('Min significance ratio maps for ROI #6');
    % sv(sigContainer);