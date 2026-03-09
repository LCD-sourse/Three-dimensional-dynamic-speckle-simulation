function [xSize, qext, qsca, albedo, g, psca0, s1tab, s2tab] = PreparMieScattering(lambda, n_water, n_sctCOMPLEX, dia, NN)
% PREPARMIESCATTERING Precompute Mie scattering parameters for Monte Carlo simulation
%
% This function calculates Mie scattering parameters including scattering amplitudes,
% efficiency factors, and anisotropy factor for given optical properties.
%
% Inputs:
%   lambda        - Wavelength in vacuum (¦Ìm)
%   n_water       - Refractive index of background medium (water)
%   n_sctCOMPLEX  - Complex refractive index of scatterer
%   dia           - Diameter of scatterer (¦Ìm)
%   NN            - Length of S1 and S2 scattering amplitude tables
%
% Outputs:
%   xSize         - Size parameter of Mie sphere
%   qext          - Extinction efficiency
%   qsca          - Scattering efficiency
%   albedo        - Single scattering albedo
%   g             - Anisotropy factor (mean cosine of scattering angle)
%   psca0         - Normalization factor for scattering phase function at forward direction
%   s1tab         - Scattering amplitude S1 table (complex values)
%   s2tab         - Scattering amplitude S2 table (complex values)

% Calculate size parameter
xSize = pi * dia * n_water / lambda;

% Calculate relative refractive index (scatterer relative to background)
mRefra = n_sctCOMPLEX / n_water;

% Call Mie scattering calculation function
% Note: The mie function should return an array with at least 9 elements
Result = mie(mRefra, xSize);

% Extract Mie scattering results
qext = Result(4);  % Extinction efficiency
qsca = Result(5);  % Scattering efficiency
qabs = Result(6);  % Absorption efficiency
qb = Result(7);    % Backscattering efficiency
g = Result(8);     % Anisotropy factor
qratio = Result(9); % Scattering efficiency ratio

% Calculate single scattering albedo
albedo = qsca / qext;

% Pre-allocate arrays for scattering amplitude tables
s1tab = zeros(1, NN);
s2tab = zeros(1, NN);
mulist = zeros(1, NN);

% Calculate scattering amplitudes for discrete scattering angles
for ii = 1:NN
    % Cosine of scattering angle (mu = cos¦È)
    % Uniform sampling from -1 to 1 (¦È from 0 to ¦Ð)
    mulist(ii) = cos(pi * (ii - 1) / (NN - 1.0));
    
    % Calculate scattering amplitudes S1 and S2 for current scattering angle
    [s1tab(ii), s2tab(ii)] = Mie_S12(mRefra, xSize, mulist(ii));
end

% Calculate normalization factor for scattering phase function at forward direction (mu=1)
% This is used for rejection sampling in the Monte Carlo simulation
psca0 = (abs(s1tab(1))^2 + abs(s2tab(1))^2) / 2;

end