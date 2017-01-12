clear all; close all; clc

%% Aufgabe 2.9 Lead-Lag Regler entwurf

Vg = 7;
a = 1/5;
b = 1
T = 1/4.2;
xi = 0.8;
t_r = 0.15;
omega_c = 1.5/t_r;
ue = 15;
phi_soll = 70 - ue;

lag_extra_phase = 10;

Gs = Vg*tf([a 1], [T^2 2*xi*T 1 0])
V_ges = 1/0.05;
V_R = V_ges / Vg;

R_1 = V_R;
L_1 = R_1*Gs;
[re_L_1 im_L_1] = nyquist(L_1, omega_c);
phi_L_1 = atan (im_L_1/re_L_1) * 180/pi;
abs_L_1 = sqrt(re_L_1^2 + im_L_1^2);

phi_dif = phi_soll - phi_L_1

% Phase um phi_dif + lag_extra_phase anheben
% Lead Regler

delta_phi = phi_dif + lag_extra_phase;
eta_lead = (tan(pi/4 - (delta_phi/2)*pi/180))^2
T_lead = 1/(sqrt(eta_lead)*omega_c)

if eta_lead < 1 
    sprintf('ja! eta lead kleiner 1')
end

R_lead = tf([T_lead 1], [T_lead*eta_lead 1])

% Lead Regler
L_2 = L_1*R_lead;
[re_L_2 im_L_2] = nyquist(L_2, omega_c);
phi_L_2 = atan (im_L_2/re_L_2) * 180/pi;
abs_L_2 = sqrt(re_L_2^2 + im_L_2^2);

phi_dif = phi_soll - phi_L_2

delta_a = 1/abs_L_2;
delta_phi = phi_dif*pi/180;

T_lag = (delta_a*sqrt(1 + tan(delta_phi)^2) - 1)/(omega_c*tan(delta_phi));
eta_lag = (omega_c*T_lag - tan(delta_phi))/(omega_c*T_lag* (1 + omega_c*T_lag*tan(delta_phi)));

R_lag = tf([T_lag 1], [T_lag*eta_lag 1])

if eta_lag > 1 
    sprintf('ja! eta lag groesser 1')
end

% Ueberpruefe
L_3 = L_2*R_lag
[re_L_3 im_L_3] = nyquist(L_3, omega_c);
phi_L_3 = atan (im_L_3/re_L_3) * 180/pi
abs_L_3 = sqrt(re_L_3^2 + im_L_3^2)

figure
line([omega_c omega_c], [25, -25])
hold on
bode(Gs, L_1, L_2, L_3)
%hold on
line([omega_c omega_c], [-90, -160])
line([5 20], [phi_soll-180,phi_soll-180])
grid on
title('2.1: Bode-Diagramm von G(s) und der offfene Regelkreise')
legend('G(s)', 'L1(s)','L2(s)','L3(s)', 'omega_c', 'phi soll');