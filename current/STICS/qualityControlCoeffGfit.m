function [coeffGfit, corrfn] = qualityControlCoeffGfit(coeffGfit, corrfn, opt)
%Quality Control. Rodrigo Migueles. Fall 2021.
% Discard correlation function fits if they have:
    cutOff = min([find(coeffGfit(:,1)<0, 1, 'first'),... % negative amplitudes or...
        find(coeffGfit(:,2)>opt.maxHalfWidth,1,'first'),... % large X width or...
        find(coeffGfit(:,3)>opt.maxHalfWidth,1,'first'),... % large Y width...
        find(coeffGfit(:,1)<coeffGfit(1,1)/4, 1, 'first'),... % or an amplitude decay beyond 1/4.
        find(abs(coeffGfit(:,1))>abs(coeffGfit(1,1).*-1.1), 1, 'first')]); % or an amplitude greater than the first.

    % Remove these fits and their corresponding correlation
    % functions:
    if cutOff >= 1
        coeffGfit = coeffGfit(1:cutOff-1,:);
        corrfn = corrfn(:,:,1:cutOff-1);
    end
end