function Omega = getAngles(u, v, degrees)
%%getAngles Calculates the angle between a vector and the horizontal.
%%Counter-clockwise from the east.
% Rodrigo Migueles Ramirez. McGill University. February 2022.

    if nargin < 3
        degrees = 0;
    end

    TanOmega = v./u;
    Theta = atan(TanOmega);

    xPyP = u>0 & v>0;
    Omega = xPyP.*Theta;

    xPyN = u>0 & v<0;
    Omega = Omega + xPyN.*(2*pi + Theta);

    xNyP = u<0 & v>0;
    Omega = Omega + xNyP.*(pi + Theta);

    xNyN = u<0 & v<0;
    Omega = Omega + xNyN.*(pi + Theta);
    
    if degrees
        Omega = rad2deg(Omega);
    end           
end