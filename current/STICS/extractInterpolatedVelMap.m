function InterpolatedVelMap = extractInterpolatedVelMap(velocityMap, opt)
InterpolatedVelMap = zeros(size(velocityMap{1,1}.iVelMap, 1), size(velocityMap{1,1}.iVelMap, 2), size(velocityMap, 2));
for m = size(velocityMap,1)
    for k = size(velocityMap, 2)
        InterpolatedVelMap(:,:,k) = velocityMap{m,k}.iVelMap;
%         imwrite(uint16(InterpolatedVelMap(:,:,k)), ...
%             [opt.path, opt.outputName, '_InterpolatedVelMap.tif'],...
%             'Compression', 'none', 'WriteMode', 'append');
    end
end