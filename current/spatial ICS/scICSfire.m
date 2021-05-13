function scICSfire(TICSScan, series, rect, ROI)
% SMV 10/2020
% Takes ICS2DScan struct array and produces a heatmap of the variable
% of choice.

    
    all_y_coords = zeros(size(TICSScan,2),1); %returns a repeating array of repeating x_coords
    all_x_coords = zeros(size(TICSScan,2),1);
    x_coords = (rect(1):1:(rect(1)+rect(3)));
    y_coords = (rect(2):1:(rect(2)+rect(4)));
    diffusions = zeros(size(TICSScan,2),1);
    partPerBA = zeros(size(TICSScan,2),1);
    BAs = zeros(size(TICSScan,2),1);
    ratio = zeros(size(TICSScan,2),1);
    
    
    for k = 1 : size(TICSScan,2)
        all_x_coords(k) = floor(TICSScan(k).location_x);
        all_y_coords(k) = floor(TICSScan(k).location_y);
        diffusions(k) = TICSScan(k).CD(1);
    end
    
    xsize = size(unique(all_x_coords),1); %returns a list of all x_coords, non repeating
    ysize = size(unique(all_y_coords),1);
    allx = zeros(ROI/2*(xsize+1),2);
    ally = zeros(ROI/2*(ysize+1),2);
    
    for x = 1:ysize
        for i = -ROI/2+1:ROI/2
            allx((x-1)*ROI+(i+ROI/2),1) = x_coords(x*ROI/2+i);
            allx((x-1)*ROI+(i+ROI/2),2) = all_x_coords(x);
        end
    end
      for x = 1:xsize
        for i = -ROI/2+1:ROI/2
            ally((x-1)*ROI+(i+ROI/2),1) = y_coords(1,x*ROI/2+i);
            ally((x-1)*ROI+(i+ROI/2),2) = all_y_coords((x-1)*xsize+1,1);
        end
    end   
    
    M = zeros(length(y_coords), length(x_coords));
    % The following loop(s) will insert each 'density' value into the appropriate position in the
    %matrix
    
    %k = ((i-1)*6)+1;
    x = 1;
    while x <= length(x_coords) %will shift to the next column as all rows are filled
        y = 1;
        xplace = find(x_coords(x) == allx(:,1));
        while y <= length(y_coords) %will allocate the y values in order
            yplace = find(y_coords(y) == ally(:,1));
            if size(xplace,1)*size(yplace,1) == 1
                index = find(allx(xplace,2) == all_x_coords(:,1) & ally(yplace,2) == all_y_coords(:,1));
                M(y,x) = diffusions(index,1);
            elseif size(xplace,1)*size(yplace,1) == 2
                nonmean = zeros(2,1);
                if size(xplace,1) == 2  
                    for i = 1:2
                        index = find((allx(xplace(i),2) == all_x_coords(:,1) & ally(yplace,2) == all_y_coords(:,1)));
                        nonmean(i,1) = diffusions(index,1);
                    end
                else
                    for i = 1:2
                        index = find((allx(xplace,2) == all_x_coords(:,1) & ally(yplace(i),2) == all_y_coords(:,1)));
                        nonmean(i,1) = diffusions(index,1);
                    end
                end
                M(y,x) = mean(nonmean, 'omitnan');
            elseif size(xplace,1)*size(yplace,1) == 4
                nonmean = zeros(4,1);
                p = 1;
                for i = 1:2
                    for j = 1:2
                        index = find((allx(xplace(i),2) == all_x_coords(:,1) & ally(yplace(j),2) == all_y_coords(:,1)));
                        nonmean(p,1) = diffusions(index,1);
                       p =  p+1;
                    end
                end
                M(y,x) = mean(nonmean, 'omitnan');
%             elseif size(xplace,1)*size(yplace,1) == 6
%                 nonmean = zeros(3,1);
%                 p = 1;
%                 if size(xplace,1) == 3
%                     for i = 1:2
%                         for j = 1:3
%                             index = find((allx(xplace(j),2) == all_x_coords(:,1) & ally(yplace(i),2) == all_y_coords(:,1)));
%                             if index
%                                 nonmean(p,1) = diffusions(index,1);
%                                 p = p+1;
%                             end
%                         end
%                     end
%                 else
%                     for i = 1:2
%                         for j = 1:3
%                             index = find((allx(xplace(i),2) == all_x_coords(:,1) & ally(yplace(j),2) == all_y_coords(:,1)));
%                             if index
%                                 nonmean(p,1) = diffusions(index,1);
%                                 p = p+1;
%                             end
%                         end
%                     end
%                 end
%                 M(y,x) = mean(nonmean);
            else
                M(y,x) = NaN;
            end
       
            %k = k + 1;
            y = y + 1;
        end
        x = x + 1;
    end
    
    figure(1)
    h = heatmap(M)
    h.GridVisible = 'off';
    h.Colormap = parula;



