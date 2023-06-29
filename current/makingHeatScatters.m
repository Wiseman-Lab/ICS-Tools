keptV = reshape(velocityMap{1,1}.goodVectors, [35, 35]);
imshow(keptV);
lvx = reshape(velocityMap{k}.vx, [size(position_x,1)*size(position_y,2),1]);
lvy = reshape(velocityMap{k}.vy, [size(position_x,1)*size(position_y,2),1]);
lv = [lvx,lvy];

lvx(isnan(lvx))=0;
lvy(isnan(lvy))=0;
outpath = 'C:\Users\malvi\OneDrive - McGill University\Data\Hayer\210108\STICSOutputs\';
outname = '210108-2-7-MYL9_reg_205px-ROI8_16x4_15x1_t100r16o16v5-2Dheatscatter.png';
heatscatter(lvx, lvy, opt.path, "HeatscatterMYL9.png");
