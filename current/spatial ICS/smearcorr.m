function unsmeared = smearcorr(imgM, timesize, transfertime)

    % By SMV 10-2020 to correct for Wiseman lab CCD image smear 
    % Mathematical derivations and algorithm for codes are available at
    % Powell, K. et al. (1999) "Restoration and frequency analysis of
    % smeared images." OSA, 38., and Knox, K. T. (2006). "Recovering
    % saturated pixels blurred by CCD image smear." 
    
    % imgM should be a background-subtracted 3D image matrix. timesize is the 'exposure time' as
    % recorded on the TIRF machine. transfertime is the amount of time that
    % the camera takes to transfer the data to its empty storage bins, and
    % should be empirically calculated. ANDOR 888 should be between 0.6 and
    % 4.33us per pixel transfer time, depending on cooldown temperature of
    % the CCD.

    % calculate alpha parameter
    timeRatio = transfertime / timesize;

    % preallocate matrix
    unsmeared = zeros(length(imgM(:,1,1)), length(imgM(1, :, 1)), length(imgM(1, 1, :)));
    
    
    for i = 1 : length(imgM(1,1,:)) %iterate through all frames of the image
      for j = 1 : length(imgM(1,:,i)) %iterate through all columns of the image
          T = sum(imgM(:, j, i)); %represents the smear contribution from the total intensity of the column
          unsmeared(:, j, i) = imgM(:, j, i) - timeRatio * T;
      end
    
end
