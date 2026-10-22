function [rho, Tair, p, c] = stdatm(z)
% STDATM  US Standard Atmosphere 1976 (troposphere, 0–11 km)
% This function is to solve for the - air density, temp of air, pressure, and
% speed of sound - at a specific given height (z). The validity of this
% model
% is from 0-11km. 
% The inputs are - rho (kg/m^3) | Tair (Kelvin) | p
% (Pascals) | c (m/s) |
% The output is z - height above sea level (m); can be vector or scalar.

% Constants (US Standard Atmosphere 1976)
T_0 = 288.15; % sea-level temp, K
p_0 = 101325; % sea-level pressure, Pa
gamma = 1.4; % ratio of specific heats
L = 0.0065; % temperature lapse rate, K/m
g0 = 9.80665; %standards defined value
R = 287.05; %J/(kg*K)
exponent = g0/(L*R);

% Error for any height above 11,000km
if any(z>11000)
    error('Model only valid up to 11km')
end
% Output
Tair = T_0 - L.*z;
p = p_0 .* ((Tair ./ T_0).^(exponent));
rho = p ./ (R.*Tair);
c = sqrt(gamma.*R.*Tair);
end
