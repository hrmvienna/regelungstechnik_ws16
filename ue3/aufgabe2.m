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

f_M = [-(R_GSM*i_GSM - u_GSM + k_GSM*w_GSM)/L_GSM;...
    w_GSM - w_P; ...
    -(d_cGSM + c_GSMP*phi_GSMP - i_GSM*k_GSM + d_vGSM*w_GSM + d_GSMP*(w_GSM - w_P))/J_GSM; ...
     -(M_ext + d_cP - c_GSMP*phi_GSMP + d_vP*w_P - d_GSMP*(w_GSM - w_P) + d_qP*w_P^2)/J_P];
 
% Parameterliste
paralist_1 = [u_GSM_r M_ext_r L_GSM    R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4e-3   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];

f_M_v = subs(f_M, paralist_1, paralist_2)

% Ruhelagen
x_R = [5.52282806081939;0.506027297604110;30.5949909202308;30.5949909202308];
i_GSM_r = x_R(1);
phi_GSMP_r = x_R(2);
w_GSM_r = x_R(3);
w_P_r = x_R(4);