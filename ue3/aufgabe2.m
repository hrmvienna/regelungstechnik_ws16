clear all; close all; clc

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

x_M = [i_GSM phi_GSMP w_GSM w_P].';
f_M = [-(R_GSM*i_GSM - u_GSM + k_GSM*w_GSM)/L_GSM;...
    w_GSM - w_P; ...
    -(d_cGSM + c_GSMP*phi_GSMP - i_GSM*k_GSM + d_vGSM*w_GSM + d_GSMP*(w_GSM - w_P))/J_GSM; ...
     -(M_ext + d_cP - c_GSMP*phi_GSMP + d_vP*w_P - d_GSMP*(w_GSM - w_P) + d_qP*w_P^2)/J_P];

% Parameterliste
paralist_1 = [u_GSM_r M_ext_r L_GSM    R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4e-3   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];
 
% Ruhelagen
x_R = [5.52282806081939;0.506027297604110;30.5949909202308;30.5949909202308];
i_GSM_r = x_R(1);
phi_GSMP_r = x_R(2);
w_GSM_r = x_R(3);
w_P_r = x_R(4); 
u_GSM_r = 5.6;
M_ext_r = 0;

% Nichtlineares System linearisieren
A2 = double(subs([ -R_GSM/L_GSM,   0,  -k_GSM/L_GSM,  0;...
  0,       0,      1,      -1;...
  k_GSM/J_GSM, -c_GSMP/J_GSM, -(d_GSMP + d_vGSM)/J_GSM, d_GSMP/J_GSM;...
  0,    c_GSMP/J_P,  d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r)/J_P], paralist_1, paralist_2))

bu2 = double(subs([1/L_GSM;0;0;0], paralist_1, paralist_2))
bd2 = double(subs([0;0;0;-1/J_P], paralist_1, paralist_2))
ct2 = double(subs([0 0 0 1], paralist_1, paralist_2))
d2 = 0;

% Abtastsystem
sys2 = ss(A2, bu2, ct2, d2)
dsys2 = c2d(sys2, Ta)
phi2 = dsys2.A;
gamma2 = dsys2.B;

% Zustandsregler für das linearisierte System berechnen

% Polvorgabe, gewuenschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
lambda0 = -4;
P = [lambda0, lambda0, lambda0, lambda0];
Pd = exp(P*Ta);

kt2 = -acker(phi2, gamma2, Pd) % Vorzeichen, Hinweis in Beispiel 8.1
% Gleichung (8.39)
g2 = 1 / ((ct2 + d2*kt2)*inv(eye(4) - phi2 - gamma2*kt2)*gamma2 + d2)


f_M_v = subs(f_M, paralist_1, paralist_2)
