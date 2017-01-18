clear all; close all; clc

%% Aufgabe 3.1.1: Zustandsreglerentwurf

% Parameter
Ta = 10e-03; % 10ms

% Lin. red. Modell aus Aufgabe 2.3
A = double([0, 1, -1; -3411/62, -54163/28520, 1/1240; ...
    6822/325, 1/3250, 54117/74750 - (7*1130385^(1/2))/7475]);
bu = double([0; 12500/713; 0]);
bd = double([0; 0; -400/13]);
ct = [0 0 1];
d = 0;
sys = ss(A,[bu, bd], ct, d)

% System abtasten
dsys = c2d(sys, Ta)

phi = dsys.A;
gamma = dsys.B(:,1); % Spalte 1

% Erreichbarkeitsmatrix
R = ctrb(phi, gamma);

if rank(R) == 3
   sprintf('System ist vollstaendig erreichbar')
else
    sprintf('System ist nicht vollstaendig erreichbar')
end

%% Aufgabe 3.1.2: diskreter Zustandsregler

syms phi_GSMP w_GSM w_P uk rk g
kt = sym('kt', [1 3]);

xk = [phi_GSMP w_GSM w_P].';
uk = kt*xk + g*rk;

dxk = phi*xk + gamma*uk;
yk = ct*xk + d*uk;

% Ueberpruefung, dxk == dxk2, Skript Kapitel 8.1
phi_g = phi + gamma*kt;
dxk2 = phi_g*xk + gamma*g*rk;

% Eigenwerte der Dynamikmatrix
eigen_phi = eig(phi);

% Polvorgabe, 1/2 wurde gewaehlt, da es kleiner als die Eigenwerte von Phi
% ist, d.h. schnellere Dynamik der gewaehlten Pole als von Phi
lambda0 = -1/2;
%P = [exp(Ta*lambda0), exp(Ta*lambda0), exp(Ta*lambda0)];
P = [-0.5, -0.5, -0.5];
kt = -acker(phi, gamma, P) % Vorzeichen, Hinweis in Beispiel 8.1
kt1 = kt(1); kt2 = kt(2); kt3 = kt(3);

% Gleichung (8.39)
g = 1 / ((ct + d*kt)*inv(eye(3) - phi - gamma*kt)*gamma + d)

phi_g = double(subs(phi_g));
% Streckenuebertragungsfunktion
sys_g = ss(phi_g, gamma*g, (ct+ d*kt), d*g);
Gz = tf(sys_g)

figure
bode(Gz, (Gz/(1+Gz)))

%% Reglerentwurf
% Anforderungen aus Aufgabe 2.5

% Sollsprung delta_r = 20 rad s^-1
delta_r = 20;

% Bleibende Reglerabweichung fuer Eingangssprung r^k = (1^k) 
e_inf_s = 0;
% Antiegszeit
t_r = 1;
% Ueberschwingen
u_e = 0;
% Stellgroessenbeschraenkung
u_gsm_min = 0;
u_gsm_max = 12;