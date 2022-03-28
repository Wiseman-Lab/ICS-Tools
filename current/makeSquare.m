function output = makeSquare(im)
%%makeSquare Converts a 2D array of different sizes into a quare one.
% RM March 2022

    if size(im,1) ~= size(im,2)
        newSize = min(size(im,1), size(im,2)) - mod(min(size(im,1), size(im,2)),2);    
        output = im(1:newSize, 1:newSize, :);
    else
        output = im;
    end
end