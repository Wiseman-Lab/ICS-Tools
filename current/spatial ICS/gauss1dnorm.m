function F = gauss1dnorm(a,data,weights);

% Used by the curve fitter to calculate values for a 1d normalized gaussian

F = (exp(-(data.^2)/(a(2)^2))).*weights;