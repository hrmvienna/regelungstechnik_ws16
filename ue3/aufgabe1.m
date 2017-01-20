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
sys = ss(A,bu, ct, d)

% System abtasten
dsys = c2d(sys, Ta)

phi = dsys.A;
gamma = dsys.B; % Spalte 1

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

% Ueberpruefung, dxk == dxk2, Skript Kapitel 8.1
phi_g = phi + gamma*kt;

% Eigenwerte der Dynamikmatrix
eigen_A = eig(A);
% eigen_A = [ -0.7225 +- 8.6575i, -0.7258]

% Polvorgabe, gewuenschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
lambda0 = -3.5;
P = [lambda0, lambda0, lambda0];
% Gewuenschte Pole des geschlossenen Kreises für das Abtastsystem
Pd = exp(P*Ta);

kt = -acker(phi, gamma, Pd) % Vorzeichen, Hinweis in Beispiel 8.1
kt1 = kt(1); kt2 = kt(2); kt3 = kt(3);

%Test
eig(phi+gamma*kt)

% Gleichung (8.39)
g = 1 / ((ct + d*kt)*inv(eye(3) - phi - gamma*kt)*gamma + d)

phi_g = double(subs(phi_g));
% Streckenuebertragungsfunktion
sys_g = ss(phi_g, gamma*g, (ct+ d*kt), d*g, Ta);
Gz = tf(sys_g)

%% Reglerentwurf
%% A: Anforderungen aus Aufgabe 2.5

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

% TODO: Anstigszeit kann durch die Polvorgabe beeinflusst werden

%% E: Sprungantwort des geschlossenen Kreises

% Geschlossener Kreis
T_ry = Gz;

figure
step(T_ry)
% Ueberschwingung einzeichnen
line([0, 7], [1, 1], 'Color', 'r')
% tr so halbwegs einzeichnen, wie Abbildung 5.2.
a = 0.76; % Wendepunkt, vom Plot abgelesen (anklicken)
line([a-t_r/2, a+t_r/2], [0, 1], 'Color','k')
line([a-t_r/2, a-t_r/2], [0, 1], 'Color','g')
line([a+t_r/2, a+t_r/2], [0, 1], 'Color','g')
title('Sprungantwort des geschlossenen Kreises L4(s)')
legend('Try', 'ue', 'tr')
grid on

stepinfo(T_ry)
