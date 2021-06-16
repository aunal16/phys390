function refractiveIndexSurfacePlotter(L, MLAT)
% This function plots the refracetive index surface plot, by default for 
% L = 4 and MLAT = 0. To get Ne and Ni, surfacedens1.m is used manually.
close all
if nargin<2
    L = 4; MLAT = 0; Ne = 7.5354e+08; Ni = 7.5351e+08; end;

q = 1.602e-19;		      % electron charge (C)
me = 9.10955854e-31;	  % electron mass (kg)
am = 1.67252e-27 / me;    % H+/e- mass ration (~1827)
eps = 8.8541878128 * 10^(-12); % vacuum permittivity

B0 = B0_calc(L, MLAT);                      % static magnetic field
% B0 = 9.7500e-07;    % bmodel.m is used manually for the input (4, 0)
wC = wC_calc(q, B0, me);                    % electron gyrofrequency
wP = wP_calc(Ne, Ni, q, eps, me, me*am);    % plasma freuquency
wH = sqrt(wP^2 + wC^2);                     % upper hybrid oscillation freq

fprintf("B0: %g\nωc: %g\nωp: %g\n", B0, wC, wP);
increment = 0.2;

wH = 2* pi * 13.5 * 1000;   % given by the instructor
theta = 0:0.01:2*pi;

figure('WindowState','maximized');
for i = 0.1 : increment : 0.5
    % plus version of refractive index
    w = wH * i;
    
    X  = wP^2 / w^2;
    YT = wC * sin(theta) / w;
    YL = wC * cos(theta) / w;
    
    n_squared_p = 1 - X ./ ( 1 - YT.^2 ./ (2*(1 - X)) + 1./(1-X) .* ...
        sqrt(YT.^4 / 4 + YL.^2 .* (1 - X).^2));
    
    mu_p_rea = real(sqrt(n_squared_p));
    mu_p_img = imag(sqrt(n_squared_p));
    mu_p_abs = abs (sqrt(n_squared_p));
    
    if i == 0.5
        plot_polar(theta, mu_p_rea, mu_p_img, mu_p_abs, increment, true, true);
    else
        plot_polar(theta, mu_p_rea, mu_p_img, mu_p_abs, increment, false,true);
    end
end

figure('WindowState','maximized');
for i = 0.1 : increment : 0.5
    % minus version of reftactive index
    w = wH * i;
    
    X  = wP^2 / w^2;
    YT = wC * sin(theta) / w;
    YL = wC * cos(theta) / w;
    
    n_squared_m = 1 - X ./ ( 1 - YT.^2 ./ (2*(1 - X)) - 1./(1-X) .* ...
        sqrt(YT.^4 / 4 + YL.^2 .* (1 - X).^2));
    
    mu_m_rea = real(sqrt(n_squared_m));
    mu_m_img = imag(sqrt(n_squared_m));
    mu_m_abs = abs (sqrt(n_squared_m));
    
    if i == 0.5
        plot_polar(theta, mu_m_rea, mu_m_img, mu_m_abs, increment, true, false);
    else
        plot_polar(theta, mu_m_rea, mu_m_img, mu_m_abs, increment, false, false);
    end
end
return

function plot_polar(theta, re, im, ab, inc, lastCheck, sign)
% This function plots real (mu), imaginary (chi) and absolute(|n|) parts
% of the refractive index (n).

rlim_end = 150;
if lastCheck
    legendLabels = cell(1, (0.5-0.1)/inc + 1);
    for i = 0.1 : inc : 0.5
        tit = strcat(sprintf("%g", i), " \omega_H");
        legendLabels{1, int8((i + 0.1) / inc)} = tit;
    end
end

subplot(3,1,1); polarplot(theta, re); rlim([0, rlim_end]); hold on;
ax = gca; ax.ThetaZeroLocation = "top"; ax.ThetaDir = "clockwise";
title("\theta vs \mu"); 
if sign
    subtitle("Real part of n^+");
else
    subtitle("Real part of n^-");
end
if lastCheck
    legend(legendLabels); end

subplot(3,1,2); polarplot(theta, im); rlim([0, rlim_end]); hold on;
ax = gca; ax.ThetaZeroLocation = "top"; ax.ThetaDir = "clockwise";
title("\theta vs \chi"); 
if sign
    subtitle("Imaginary part of n^+");
else
    subtitle("Imaginary part of n^-");
end
if lastCheck
    legend(legendLabels); end

subplot(3,1,3); polarplot(theta, ab); rlim([0, rlim_end]); hold on;
ax = gca; ax.ThetaZeroLocation = "top"; ax.ThetaDir = "clockwise";
title("\theta vs |n|");
if sign
    subtitle("Magnitude of n^+");
else
    subtitle("Magnitude of n^-");
end
if lastCheck
    legend(legendLabels); end
return

function B0 = B0_calc(L, MLAT)
% lambda is geomagnetic latitude (MLAT) in degrees
% L is L-shell
lambda = MLAT;
rat = 1 / (L * (cosd(lambda))^2);       % Ratio of R_0/R
B0 = 0.312 * 10^(-4) * (rat)^3 * sqrt(1 + 3*(sind(lambda))^2);
return

function wC = wC_calc(q, B, m)
% This function calculates particle gyrofrequency (a.k.a cyclotfon freq.)

% q is particle charge  [Coulomb]
% B0 is B-field strength [Wb/m^2]
% m is particle mass [Kg]
wC = q * B / m;
return

function wP = wP_calc(Ne, Ni, q, eps, me, mi)
% This function calculates thec  complete particle plasma frequency,
% wp = wpi + wpe
wPe_squared = Ne * q^2 / (eps * me);
wPi_squared = Ni * q^2 / (eps * mi);

wP = sqrt(wPe_squared + wPi_squared);
return

