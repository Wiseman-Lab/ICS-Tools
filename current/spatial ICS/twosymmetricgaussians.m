function F = twosymmetricgaussians(a,data,weights);

% Used by the curve fitter to calculate values for a 2d gaussian

X = data(:,1:size(data,2)/2);
Y = data(:,size(data,2)/2+1:end);


F = ((a(1)*exp(    -((X-a(2)).^2/(a(7)^2)+(Y-a(3)).^2/(a(7)^2))   ) + a(8))  + (a(4)*exp(    -((X-a(5)).^2/(a(7)^2)+(Y-a(6)).^2/(a(7)^2))) )) .*weights;

