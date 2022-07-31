function [c_mask]=circle(ix,iy,cx,cy,r)

% By M. Bach
% Draws a circle with centre cx,cy in image ix,iy with radius r
% [c_mask]=circle(ix,iy,cx,cy,r)

[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
c_mask=((x.^2+y.^2)<=r^2);