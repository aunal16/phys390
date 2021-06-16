function dispersionRelationPlotter(L, MLAT, N0)

close all;

if nargin<1
    L = 1;
    MLAT = 20;
    N0 = 2.1673e+14;
    N0_p = 1.5611e+05;
end

R_earth = 6371*1000;
q = 1.602e-19;		      % electron charge (C)
me = 9.10955854e-31;	  % electron mass (kg)
am = 1.67252e-27 / me;    % H+/e- mass ration (~1827)
eps = 8.8541878128 * 10^(-12); % vacuum permittivity
c = 3*10^8;

height = L * R_earth;
B0 = B0_calculate(L, MLAT);
wc = wCyc_calculate(q, B0, me); % electron gyrofreq.
wpe = wpe_calculate(N0, q, me, eps); % plasma oscillation (electron) freq.
wpi = wpe_calculate(N0_p, q, am*me, eps); % plasma oscillation (electron) freq.
wp = sqrt(wpe^2+wpi^2);
% wh = sqrt(wp^2+wc^2); % upper hybrid frequency
wh = sqrt(wpe^2+wc^2); % upper hybrid frequency

fprintf("At %g km and %g ° MLAT\n\tB0 is %g\n\twc is %g\n\twpe is %g\n\twpi is %g\n\twp is %g\n\twh is %g\n", height/1000, MLAT, B0, wc, wpe, wpi, wp, wh);

%% theta = 0
omega1 = linspace(0, wc-100, 500);
wR = sqrt(wpe^2+wc^2/4)+wc/2;
wL = sqrt(wpe^2+wc^2/4)-wc/2;
omega2 = linspace(wR, wpe*4, 500);
omega3 = linspace(wL, wpe*4, 500);

k1 = omega1/c.*sqrt(1-wpe^2./(omega1.*(omega1-wc)));
k2 = omega2/c.*sqrt(1-wpe^2./(omega2.*(omega2-wc)));
k3 = omega3/c.*sqrt(1-wpe^2./(omega3.*(omega3+wc)));
k4 = linspace(0,10);

subplot(3,1,1)
plot(k3, omega3, "b")
hold on;
subplot(3,1,1)
plot(k2, omega2, "r")
ylim([0.81 0.9]*10^9)
xlim([0,3])
hold on;
subplot(3,1,1)
yline(wpe, "k","\omega_{p}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "middle");
subplot(3,1,1)
yline(wR, "--","\omega_{R}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "top");
subplot(3,1,1)
yline(wL, "--","\omega_{L}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "bottom");
hold on;

grid on;
ylabel("\omega");
title("Dispersion relation for \theta=0");
legend("LH","RH","\omega_p")

subplot(3,1,2)
plot(k3, omega3,"b")
hold on;
plot(k2, omega2,"r")
xlim([0,10])
hold on;
plot(k4, k4*c)
yline(wpe, "k","\omega_{p}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "middle");
yline(wc, "m","\omega_{c}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "top");
% hold on;
% yline(wR, "--","\omega_{R}", "LabelHorizontalAlignment", "Center", ...
%     "LabelVerticalAlignment", "top");
% hold on;
% yline(wL, "--","\omega_{L}", "LabelHorizontalAlignment", "Center", ...
%     "LabelVerticalAlignment", "bottom");
% hold on;

grid on;
ylabel("\omega");
legend("LH","RH","\mu=1","\omega_p", "\omega_c","Location", "Northwest")

subplot(3,1,3)
plot(k1, omega1, "r"); hold on;
yline(wc,"m", "\omega_{c}", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "middle");
xlim([0,10])

grid on;
xlabel("k");
ylabel("\omega");
legend("RH","\omega_c","Location", "Southeast")

%% theta = 90
figure
omega4 = linspace(wL, wh-10, 2000);
omega5 = linspace(wpe, wpe*4, 2000);
omega6 = linspace(wR, wpe*4, 2000);

k4 = omega4./c.*sqrt(1-wpe^2./omega4.^2.*((omega4.^2-wpe^2)./(omega4.^2-wpe^2-wc^2)));
k5 = omega5./c.*sqrt(1-wpe^2./omega5.^2);
k6 = omega6./c.*sqrt(1-wpe^2./omega6.^2.*((omega6.^2-wpe^2)./(omega6.^2-wpe^2-wc^2)));
k7 = linspace(0,10);

subplot(2,1,1)
plot(k5, omega5); hold on;
plot(k6, omega6, "r"); hold on;
plot(k4, omega4, "r"); hold on;
xlim([0 1.5]);
ylim([8.1 8.9]*10^8);
yline(wR, "--", "\omega_R", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "top");
yline(wh, "--", "\omega_H", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "middle");
yline(wL, "--", "\omega_L", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "bottom");

ylabel("\omega")
title("Dispersion relation for \theta=\pi/2");
grid on;
legend("O","E","E");

subplot(2,1,2)
plot(k5, omega5); hold on;
plot(k6, omega6,"r"); hold on;
plot(k4, omega4, "r"); hold on;
plot(k7, k7*c);
xlim([0 10]);
yline(wh, "--", "\omega_H", "LabelHorizontalAlignment", "Center", ...
    "LabelVerticalAlignment", "top");

ylabel("\omega")
xlabel("k");
grid on;
legend("O", "E", "E", "\mu=1", "Location", "Northwest");
return

function B0 = B0_calculate(L, lambda)
% lambda is geomagnetic latitude (MLAT) in degrees
% L is L-shell
rat = 1 / (L * (cosd(lambda))^2);       % Ratio of R_0/R
B0 = 0.312 * 10^(-4) * (rat)^3 * sqrt(1 + 3*(sind(lambda))^2);
return

function cycFreq = wCyc_calculate(q, B0, m)
% This function calculates particle gyrofrequency (a.k.a cyclotfon freq.)

% q is particle charge  [Coulomb]
% B0 is B-field strength [Wb/m^2]
% m is particle mass [Kg]

cycFreq = q * B0 / m;
return

function wpe = wpe_calculate(N0, q, me, eps)
% This function calculates electron plasma oscillation.

wpe = sqrt(N0*q^2/(me*eps));
return