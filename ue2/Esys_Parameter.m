% Parameterfiles fuer das Simulink Modell des elektrischen Systems.
close all;
clear all;
clc;

% Widerstaende R1 und R2
R_1 = 625;
R_2 = 3500;

% Widerstandsverhaeltnis R4
K = 3;

% Kapazitaeten C1 und C2
C1_ref  = 1e-6;
UC1_ref = -10;
kC1     = 800e-9;
C2      = 1e-6;

% Festlegung der Ruhelage
UeR = 5;
UsR = 5;

%% Eingangsspannung und Stoerspannung
% Sprungfoermige Eingangsspannung
te = .05;      % Einschaltzeitpunkt
Ue = 7;      % Endwert des Sprunges

% Sinusfoermige Eingangsspannung
Usinus = 1;    % Amplitude des Sinus    
we     = 1000; % Winkelfrequenz des Sinus   

% Sprungfoermige Stoerspannung
ts = .1;     % Einschaltzeitpunkt   
Us = 4;     % Endwert des Sprunges
% Rampenförmige Störspannung
tsr = 0.15; %Einschaltzeitpunkt
dr = 10;    % Steigung


%% Modell des linearisierten Systems

R1 = R_1;
R2 = R_2;

% Aufgabe 1.4.4
% -------------
% Errechnen Sie die Ruhelage des Systems!
UC1R = -(K*R2*UeR - R1*UsR - 2*R2*UeR + K*R1*UsR + K*R2*UsR)/(R1 + 2*R2);
UC2R = (R2*UeR + R1*UsR + R2*UsR)/(R1 + 2*R2);
UaR  = K*UC2R;

% Errechnen Sie fuer die gegebene Ruhelage die Ableitungen der in C1 
% gespeicherter Ladung Q1 nach der anliegenden Spannung!
dQ1  = C1_ref + UC1R*kC1 - UC1_ref*kC1; % Erste Ableitung von Q1 nach UC1
ddQ1 = kC1; % Zweite Ableitung von Q1 nach UC1
C1 = dQ1;

% Ergaenzen Sie die Systemmatrix A (Asys), die Eingangsvektoren bu (busys)
% und bd (bdsys) sowie den Ausgangsvektor c (csys) und den Durchgriff d (dsys)
% fuer das linearisierte System!

A11 = -(R1 + R2)/(C1*R1*R2);
A12 = -(K*R1 - R1 + K*R2)/(C1*R1*R2);
A21 = 1/(C2*R2);
A22 = (K - 2)/(C2*R2);

A  = [A11,A12;A21,A22];
bu = [1/(C1*R1);0];
bd = [0;1/(C2*R2)];
c  = [0, K];
du = [0];
dd = [0];

%% Regler

R_PI = tf([{[0.299035517406156 369.804412414865]}], [1 0]);
R_PID = tf([{[0.00925481711944263 9.98103353193355 726.216412490922]}], [{[0.0261814020857752 1 0]}]);

