clear all; close all; clc

%% Aufgabe 3.3: Zeitdiskreter PI-Zustandsregler

% Lin. red. Modell aus Aufgabe 2.3
A = double([0, 1, -1; -3411/62, -54163/28520, 1/1240; ...
    6822/325, 1/3250, 54117/74750 - (7*1130385^(1/2))/7475]);
bu = double([0; 12500/713; 0]);
bd = double([0; 0; -400/13]);
ct = [0 0 1]; d = 0;
sys = ss(A,bu, ct, d);

% Abtastsystem
Ta = 10e-03; % 10ms
dsys = c2d(sys, Ta);
phi = dsys.A;
gamma = dsys.B;

%% PI-Zustandsregler - Skript Kapitel 8.2, (8.52)
AI = [phi,[0,0,0]'; -ct, 1];
bI = [gamma; 0];

% Gewuenschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
lambda0 = -4; PI = [lambda0, lambda0, lambda0, lambda0];
PdI = exp(PI*Ta); % gew. Pole des geschl. Kreises für das Abtastsystem

% Polvorgabe
ke = -acker(AI,bI,PdI);

%Test
eig(AI+bI*ke)

%Berechnung der Regleranteile gemaess Skriptum
kI = ke(4);
kP = 1/(ct*inv(eye(3)-phi)*gamma);
kx = kP*ct+ke(1:3);
