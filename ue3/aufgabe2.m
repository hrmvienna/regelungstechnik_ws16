clear all; close all; clc

%% Theoriefrage

% TODO: Erfüllt der Regler die Anforderungen?

% Warum hat der Regler beim nichtlinearen System eine andere
% Reglerabweichung als beim linearen System? 

% Auch bei der Simulation 2.6 erkennt man eine bleibende Regelabweichung
% fuer die beiden linearisierten Regler. 

% Der lineare Regler verhält sich auf einer nicht-linearen Strecke
% nicht ideal und es bleibt eine Regelabweichung. warumm?

%% Aufgabe 3.2: Implementierung des Zustandsreglers aus Aufgabe 3.1

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

g = 0.165707650406087;
kt = [1.894974536721893,-0.520134683695267,0.478777536929102];

%% Zustandsregler am nichtlinearen Model aus Aufgabe 2.2

% Zustandsgroessen
syms x_m x_M i_GSM phi_GSM w_GSM phi_P w_P phi_GSMP M_ext u_GSM_r M_ext_r
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
syms J_P d_cP d_vP d_qP
syms c_GSMP d_GSMP phi_GSMP

% Ruhelagen
x_R = [5.52282806081939;0.506027297604110;30.5949909202308;30.5949909202308];
i_GSM_r = x_R(1);
phi_GSMP_r = x_R(2);
w_GSM_r = x_R(3);
w_P_r = x_R(4); 
u_GSM_r = 5.6;
M_ext_r = 0;

%% Parameter fuer NLIN - System aus Aufgabe4

u_GSM_r = 5.6;
M_ext_r = 0;

% Parameter
L_GSM   = 1.4e-3;
R_GSM   = 0.46;
k_GSM   = 0.1;
J_GSM   = 12.4e-03;
d_cGSM  = 0.152;
d_vGSM  = 1.8e-03;
J_P     = 32.5e-03;
d_cP    = 0.169;
d_vP    = 2.7e-03;
d_qP    = 1e-04;
c_GSMP  = 0.6822;
d_GSMP  = 1e-05;

% Anfangsbedingungen fuers Nicht-Lineare System
i_GSM_0 = i_GSM_r;
phi_GSMP_0 = phi_GSMP_r;
w_GSM_0 = w_GSM_r;
w_P_0 = w_P_r;
x_r = [i_GSM_r; phi_GSMP_r; w_GSM_r; w_P_r];