function [ rho_x ] = CalculateDensity( pressure,temperature )
%UNTITLED Summary of this function goes here
%   Calculates the density rho (vgl. paper yue)
%   Input-Druck in Pa
%   Input-Temperatur in Kelvin
%   Input in Vektorform

M = 28.949; % Molare Masse in g/mol
R = 8.3143; % Molare/Allgemeine Gaskonstante [J/mol*K]
rho_x = zeros(size(pressure));

% translation of the units (the following equations only works with Celsius & Mpa)
p = pressure/(10^6); % Pa -> Mpa
t = temperature-273.15; % K -> °C

for i=1:size(p)

    % berechnet Z mit P[MPa] und T[°C]
    Z(i) = 1.00001 - 5.8057*10^(-3)*p(i) + 2.6402*10^(-4)*p(i)^2 - 3.3297*10^(-7)*t(i) + 1.2420*10^(-4)*p(i)*t(i) - 2.0158*10^(-6)*p(i)^2*t(i) + 2.4925*10^(-9)*t(i)^2 - 6.2873*10^(-7)*p(i)*t(i)^2 + 5.4174*10^(-9)*p(i)^2*t(i)^2;

    % berechnet Dichte in g/m^3
    rho_x(i) = ((p(i)*10^6*M)/(R*(t(i)+273.15)*Z(i)));

    % berechnet Dichte in kg/m^3
    rho_x(i) = rho_x(i)/1000; 
end
end

