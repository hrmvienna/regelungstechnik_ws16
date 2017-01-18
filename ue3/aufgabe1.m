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
R = [gamma, phi*gamma, phi*phi*gamma];

if rank(R) == 3
   sprintf('System ist vollstaendig erreichbar')
else
    sprintf('System ist nicht vollstaendig erreichbar')
end

%% Aufgabe 3.1.2: diskreter Zustandsregler

syms delta_phi_GSMP delta_w_GSM delta_w_P delta_uk delta_rk

delta_x_red = [delta_phi_GSMP delta_w_GSM delta_w_P].';

% Polvorgabe
P = [1, 2, 3];
kt = acker(phi, gamma, P)