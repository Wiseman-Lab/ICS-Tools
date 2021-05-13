function F = twogauss2dwxyassymetric(a,data,weights);

% Used by the curve fitter to calculate values for a 2d gaussian

X = data(:,1:size(data,2)/2);
Y = data(:,size(data,2)/2+1:end);
Xprim=X.*cos(a(7))-Y.*sin(a(7));
Yprim=X.*sin(a(7))+Y.*cos(a(7));

F = ((a(1)*exp(    -((Xprim-a(5)).^2/(a(2)^2)+(Yprim-a(6)).^2/(a(3)^2))   ) + a(4))  + (a(8)*exp(    -((X-a(11)).^2/(a(9)^2)+(Y-a(12)).^2/(a(10)^2))) )) .*weights;

