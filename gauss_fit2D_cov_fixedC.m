function  [parFit, ssr] = gauss_fit2D_cov_fixedC(img,coordinates, center0)

%
% [parFit, ssr] = gauss_fit2D_cov_fixedC(img,coordinates, center0)
%____________________________________________
% Fit a two-dimensional gaussian to a subimage.
%
% INPUT
%
% img is the subimage to be fit 
%
% coordinates is a cell array with the coordinates of the pixels of the
% image:
% coordinates{1}: x coordinates
% coordinates{2}: y coordinates
%
% center0 are the initial guesses for the center of the gaussian
% center0(1) = x-coordinate
% center0(2) = y-coordinate
%
% OUTPUT
%
% parFit are the output of the fit and in particular:
% parFit(1) = amplitude of gaussian
% parFit(2) = sigma of gaussian
% parFit(3) = Bkg
% parFit(4) = xcoord of gaussian center
% parFit(5) = ycoord of gaussian center
%
% SSR is the squared sum of the residuals between the image and the fit.
%____________________________________________
%
% 2014. Davide Mazza. San Raffaele Scientific Institute. Milan. Italy
%____________________________________________



x = double(coordinates{1});
y = double(coordinates{2});

[px,py] = meshgrid(x,y);

% Starting values


par0(1) = max(img(:)) -  mean(img(:)); % max intensity
par0(2) = length(x)/8; % sigma
par0(3) = 0; % background
par0(4) = center0(1);
par0(5) = center0(2);


% boundaries for the fit
LB = [0,0,0, center0(1) - 2, center0(2) -2];
UB = [65536,4,65536, center0(1) + 2, center0(2) + 2];

% fit options
options = optimset('lsqnonlin');
options.Display = 'none';

% fit
parFit = lsqnonlin(@(P)objfun(P,px,py,img),par0,LB,UB,options);
residuals = objfun(parFit,px,py,img);
ssr = sum(residuals.^2);


% Object function
% --------------------
function residuals = objfun (par, x, y, img, center)
sigma_matrix = [1/par(2), 0; 0, 1/par(2)];
exponent = [x(:) - par(4), y(:) - par(5)]*sigma_matrix;
model = par(3) + par(1)*exp(-sum(exponent.*exponent,2)/2);
residuals = model - img(:);



 





