%% Generating correlation panels to visualize the significance 
series = readFileToStack(opt.filePath{1});
% Run stics_vectormapping to generate corrfn, then open stics.m and
% run the following:
imgser = regionanalyse;
upperTauLimit = min(opt.tauLimit, size(regionanalyse,3));

% Then run parts of stics.m to generate timecorr WITHOUT checking for
% significance. This way timecorr will be bigger than corrfn. Then run the
% following:

figure
hold on
for l = 1:opt.TOIsize
    subplot(3, opt.TOIsize, l); imagesc(regionanalyse(:,:,l)); axis image; hold on
    subplot(3, opt.TOIsize, l+opt.TOIsize); imagesc(timecorr(:,:,l)); axis image; hold on
end

for n = 1:size(corrfn, 3)
    subplot(3, opt.TOIsize, 2*opt.TOIsize+n); imagesc(corrfn(:,:,n)); axis image; hold on
end
hold off;
