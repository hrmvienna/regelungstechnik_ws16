% Parameterfiles fuer das Simulink Modell des elektrischen Systems 
% (Aufgabe 1.4 und 1.5)
close all;
clear all;
clc;

% Parameter
R1      = 625;    % Widerstand in Ohm
R2      = 3500;   % Widerstand in Ohm
K       = 3;      % Spannungsverstaerkung
C1_ref  = 1e-6;   % Referenz-Kapazitaetswert auf Kennlinie von C1 in F
UC1_ref = -10;    % Referenz-Spannungswert auf Kennlinie von C1 in V
kC1     = 800e-9; % Steigung der Kennlinie von C1 in F/V
C2      = 1e-6;   % Kapazitaet in F

% Festlegung der Ruhelage in Volt
UeR = 5;
UsR = 5;

%% Eingangsspannung und Stoerspannung (Aufgabe 1.5)
% Sprungfoermige Eingangsspannung
te = 1;      % Einschaltzeitpunkt in Sekunden
Ue = 7;      % Endwert des Sprunges in Volt

% Sinusfoermige Eingangsspannung
Usinus = 0.1; % Amplitude des Sinus in Volt
we     = 10;  % Winkelfrequenz des Sinus in rad/s

% Sprungfoermige Stoerspannung
ts = 3;     % Einschaltzeitpunkt in Sekunden
Us = 1;     % Endwert des Sprunges in Volt

%% Modell des linearisierten Systems

% Aufgabe 1.4.4
% -------------
% Errechnen Sie die Ruhelage des Systems!
UC1R = -(K*R2*UeR - R1*UsR - 2*R2*UeR + K*R1*UsR + K*R2*UsR)/(R1 + 2*R2);
UC2R = (R2*UeR + R1*UsR + R2*UsR)/(R1 + 2*R2);
UaR  = K*UC2R;

% Errechnen Sie fuer die gegebene Ruhelage die Ableitungen der in C1
% gespeicherter Ladung Q1 nach der anliegenden Spannung!
dQ1  = C1_ref + UC1*kC1 - UC1_ref*kC1; % Erste Ableitung von Q1 nach UC1
ddQ1 = kC1; % Zweite Ableitung von Q1 nach UC1

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
%% Uebertragungsfunktion G (Eingang u -> Ausgang y) und Uebertragungsfunktion Gd (Stoerung d -> Ausgang y)

% Aufgabe 1.4.4
% -------------
% Bestimmen Sie zunaechst eine MISO- oder zwei SISO-Zustandsraum-
% darstellungen mittels ss(). Anschliessend koennen Sie in beiden Faellen
% die gesuchten Uebertragungsfunktionen (G und Gd) mittels tf() bestimmen.
syms s
G  = (1.371e06) / (s^2 + 1600*s + 9.959e05);
Gd = (857.1 s + 1.616e06) / (s^2 + 1600 s + 9.959e05);

%% Verstaerkungsfaktor V, Daempfungsgrad xi und Zeitkonstante T von G

% Aufgabe 1.4.5
% -------------
% Bestimmen Sie den Verstaerkungsfaktor V, den Daempfungsgrad xi und die
% Zeitkonstante T der Uebertragungsfunktion G.

% V  = 1.371e06;
% T  = 0.0010;
% xi = 0.8016;

%% Bodediagramme der Uebertragungsfunktionen G und Gd

% Aufgabe 1.4.6
% -------------
% Zeichnen und interpretieren Sie die Bodediagramme der beiden
% Uebertragungsfunktionen G und Gd. Verwenden Sie dazu den Befehl bode().