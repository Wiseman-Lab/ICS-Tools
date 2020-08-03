function [significant] = correlationSignificance(corr);

%Decides if the maximum of a STICS correlation function is significant
%using the criteria in Ji and Danuser, J. Microsc., 2005

allMaxima = imregionalmax(corr);
threshold = 0.5;

% %if the max isn't much higher than the mean, forget it
% if (max(corr(:))/mean(corr(:))) < 1.2
%     significant = 0;
%     return
% end

% If there's only one maximum, it's significant
% otherwise, test the global max against the others
if sum(allMaxima(:))==1
    significant = 1;
else
    [corrMax corrMaxLocation] = max(corr(:));
    [corrMaxLocationi corrMaxLocationj] = ind2sub(size(corr),corrMaxLocation);

    if mod(size(corr,1),2)==0
        corr00 = size(corr,1)/2+1;
    else
        corr00 = round(size(corr,1)/2);
    end

    radius = max(5,0.5*sqrt( (corrMaxLocationi-corr00)^2 + (corrMaxLocationj-corr00)^2 ));
    mask = circle(size(corr,1),size(corr,2),corrMaxLocationi,corrMaxLocationj,radius);
    referenceScore = mean(mean(corr(mask)));
    significanceValues = corr(allMaxima) - referenceScore;
    significanceRatios = significanceValues/max(significanceValues);
    significanceRatios = sort(significanceRatios);

    % The 3rd last ratio in the list is the 3rd highest maximum... if it
    % "passes" the test, then all of the others will as well, and we
    % discard the correlation function.  This is for functions where two local maxima reflect the same correlation peak:
    %                  /\
    %                 /  \/\
    %                /      \  This should count as a valid function to fit!
    %               /        \
    %              /          \
    %             /            \
    % ___________/              \_____________________
    %If there are only 2 maxima, then compare the global with the other one
    if length(significanceRatios)==2
        if significanceRatios(end-1) < threshold
            significant = 1;
        else
            significant = 0;
        end
    else
        if significanceRatios(end-2) < threshold
            significant = 1;
        else
            significant = 0;
        end
    end

end