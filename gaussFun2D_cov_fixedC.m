function FitImg = gaussFun2D_cov_fixedC(par, coordinates)

% FitImg = gaussFun2D_cov_fixedC(par, coordinates)
%____________________________________________
% Compute a two dimensional gaussian
% 
% INPUT
% Coordinates cell array: 
% Coordinates{1}: x coordinates
% Coordinates{2}: y coordinates
%
% par(1) = amplitude of gaussian
% par(2) = sigma of gaussian
% par(3) = Bkg
% par(4) = xcoord of gaussian center
% par(5) = ycoord of gaussian center
%
% OUTPUT
% FitImg is the resulting image
%____________________________________________
%
% 2014. Davide Mazza. San Raffaele Scientific Institute. Milan. Italy
%____________________________________________


x = double(coordinates{1});
y = double(coordinates{2});
[px,py] = meshgrid(x,y);


exponent = ((px - par(4))/par(2)).^2 + ((py - par(5))./par(2)).^2;
FitImg = par(3) + par(1).*exp(-exponent*0.5);
end
