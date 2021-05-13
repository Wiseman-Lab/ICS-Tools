function F = gauss2dPSFremoved(a,data,weights);

% Used by the curve fitter to calculate values for a 2d gaussian

X = data(:,1:size(data,2)/2);
Y = data(:,size(data,2)/2+1:end);

F = ((a(1)*(exp(-(X.^2+Y.^2)/(a(2)+a(3))) + exp((X.^2+Y.^2)/a(3))))+a(4)).*weights;