function velocityMap = timeFilter(velocityMap,TOI_vector,n_sigma,varargin)

%This function filters false fit vectors that consist for only a single
%time frame
%   There are two methods of filtering that are applied:
%
%   The first method filters out vectors if the ones preceding them and
%   following them are both NaNs - these vectors themselves are set to NaN
%
%   The second method filters out vectors if the difference between that
%   vector and the one preceding it exceeds a certain threshold set by the
%   user in opt.sigmaTimeFilter. A lower threshold should be applied for
%   noisier data


%load([opt.path 'VelocityMap' opt.outputName '.mat']); %in case opt not in
%workspace... should have opt as an input variable.


for k = 1:length(TOI_vector)
    
    %Calling up vectors in x and y from first and second TOI, and setting bad
    %vectors to NaN
    
    TOI_pair = [TOI_vector(k), TOI_vector(k)+1]; %selects a pair of subsequent TOIs to analyze
    m=size(velocityMap{1}.vx,1);
    n=size(velocityMap{1}.vx,2);
    
    vx_1 = velocityMap{TOI_pair(1)}.vx(:);
    vy_1 = velocityMap{TOI_pair(1)}.vy(:);
    bad_1 = find(velocityMap{TOI_pair(1)}.goodVectors == 0);
    vx_1(bad_1) = NaN;
    vy_1(bad_1) = NaN;
    
    vx_2 = velocityMap{TOI_pair(2)}.vx(:);
    vy_2 = velocityMap{TOI_pair(2)}.vy(:);
    bad_2 = find(velocityMap{TOI_pair(2)}.goodVectors == 0);
    vx_2(bad_2) = NaN;
    vy_2(bad_2) = NaN;
    
    %Gives linear indices exclusively for entries which are not NaNs in both the
    %first and the second TOI
    index = zeros(size(vx_1));
    for i=1:length(vx_1)
        if ~isnan(vx_1(i)) && ~isnan(vx_2(i));
            index(i)=1;
        end
    end
    index = find(index);
    
    %Take the magnitude of the differences in vx and vy between the two
    %subsequent frames
    diffx = abs(vx_1(index)-vx_2(index));
    diffy = abs(vy_1(index)-vy_2(index));
    
    diff = sqrt(diffx.^2+diffy.^2);
    
    
    %Plot histogram of the magnitude of the difference of vectors in two subsequent TOIs
    %figure()
    %hist = histogram(diff);
    %title(['Difference in velocity magnitudes of TOI ',num2str(TOI_pair(1)),' and ',num2str(TOI_pair(2))],...
    %    'interpreter','latex','fontsize',16)
    
    for i = 1:length(index)
        if diff(i) > (mean(diff)+ n_sigma*(std(diff)))
            if k == 1
                vx_1(index(i)) = NaN;
                vy_1(index(i)) = NaN;
                
                vx_2(index(i)) = NaN;
                vy_2(index(i)) = NaN;
            else
                vx_2(index(i)) = NaN;
                vy_2(index(i)) = NaN;
            end
        end
    end
    
    %Remove vectors that are both preceded and followed by NaNs
    if k > 1
        TOI_before = TOI_vector(k)-1;
        vx_0 = velocityMap{TOI_before}.vx(:);
        
        for j = 1:length(vx_1)
            if isnan(vx_0(j)) && isnan(vx_2(j))
                vx_1(j) = NaN;
                vy_1(j) = NaN;
            end
        end
    end
    
    %Save the filtered vector maps into the velocityMap cell
    velocityMap{TOI_pair(1)}.vx = reshape(vx_1, m, n);
    velocityMap{TOI_pair(1)}.vy = reshape(vy_1, m, n);
    velocityMap{TOI_pair(2)}.vx = reshape(vx_2, m, n);
    velocityMap{TOI_pair(2)}.vy = reshape(vy_2, m, n);
    
end
end
