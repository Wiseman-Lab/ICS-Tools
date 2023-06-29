function F = gauss2dellipticalsimp(a,data,weights);

% Used by the curve fitter to calculate values for a 2d gaussian

X = data(:,1:size(data,2)/2);
Y = data(:,size(data,2)/2+1:end);

%F = (a(1)*exp(    -((X-a(5)).^2/(a(2)^2)+(Y-a(6)).^2/(a(3)^2))   ) + a(4)).*weights;
%F = (a(1)*exp(    -(1/2)*((X-a(5)).^2*(cos(a(7))^2/(a(2)^2)+sin(a(7))^2/(a(3)^2)  +  (Y-a(6)).^2*(sin(a(7))^2/(a(2)^2)+cos(a(7))^2/(a(3)^2))  +   (Y-a(6))*(X-a(5))*(sin(2*a(7))/(a(3)^2)-sin(2*a(7))/(a(2)^2)))))+ a(4)).*weights;
F = ((a(1)*exp(    -(1/2)*((X-a(5)).^2*a(2)  +  (Y-a(6)).^2*a(3)  +   (Y-a(6))*(X-a(5))*a(7))))+ a(4)).*weights;