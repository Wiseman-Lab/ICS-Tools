function plotACFSurface(corrfn, pixel_size, CoeffGfit)
%PLOTACFSURFACE Plots the surface of the correlation function and the
%fitted Gaussian, if provided.
%   Rodrigo Migueles Ramirez, McGill University, Feb 2022


[meshX, meshY] = meshgrid(-((size(corrfn,2)-1)/2)*pixel_size:pixel_size:((size(corrfn,2)-1)/2)*pixel_size,-((size(corrfn,1)-1)/2)*pixel_size:pixel_size:((size(corrfn,1)-1)/2)*pixel_size);

A = CoeffGfit(1);
Sx = CoeffGfit(2);
Sy = CoeffGfit(3);
Offset = CoeffGfit(4);
X_0 = CoeffGfit(5);
Y_0 = CoeffGfit(6);
Z = A.*exp(-(((meshX-X_0).^2)/(2*Sx.^2) + ((meshY-Y_0).^2)/(2*Sy.^2))) + Offset;

s = surf(ICF3D, meshX, meshY, corrfn(:,:,tau+1)); axis tight
s.EdgeColor = 'none';
s.FaceAlpha = 0.75;

shading interp
hold(ICF3D, 'on');
if CoeffGfit(1) > 0
    m = mesh(ICF3D, meshX, meshY, Z);
    m.EdgeColor = [0 0 0];
    m.EdgeAlpha = 0.5;
    m.FaceColor = [1,1,1];
    m.FaceAlpha = 0.25;

    [~, c] = contour(ICF3D, meshX, meshY, corrfn(:,:,tau+1));
end