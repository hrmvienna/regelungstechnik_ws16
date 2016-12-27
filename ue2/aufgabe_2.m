clear all; close all; clc

%% Aufgabe 2.2.1: Zustandstransformation

% Zustandsgroessen
syms x_m x_M i_GSM phi_GSM w_GSM phi_P w_P phi_GSMP M_ext
% GSM
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
% P
syms J_P d_cP d_vP d_qP
% GSMP
syms c_GSMP d_GSMP phi_GSMP

% Zustandsvektoren
x_m = [i_GSM phi_GSM w_GSM phi_P w_P].';
x_M = [i_GSM phi_GSMP w_GSM w_P].';

u_m = u_GSM;
d_m = M_ext;

%phi_GSMP = (phi_GSM - phi_P);

M_P = d_cP + d_vP * w_P + d_qP * (w_P)^2 + M_ext;
M_rGSM = d_cGSM + d_vGSM * w_GSM;
M_kopp = (w_GSM - w_P) * d_GSMP + (phi_GSM - phi_P) * c_GSMP;
M_GSM = k_GSM * i_GSM

% Systemgleichungen fuer System _m, praefix d_ 
d_i_GSM = (u_GSM - R_GSM * i_GSM - k_GSM * w_GSM)/L_GSM;
d_phi_GSM = w_GSM;
d_w_GSM = (M_GSM - M_rGSM - M_kopp)/J_GSM;
d_phi_P = w_P;
d_w_P = (M_kopp - M_P)/J_P;

f_m = [d_i_GSM d_phi_GSM d_w_GSM d_phi_P d_w_P].';
d_x_m = f_m

% Zustandstransformation
T =  [1 0 0 0 0;    0 1 0 -1 0; 0 0 1 0 0;  0 0 0 0 1];

f_M = T*f_m;
d_x_M = f_M

u_M = u_m;
d_M = d_m;

%% Aufgabe 2.2.2:  Linearisieren um die Ruhelage

syms u_GSM_r M_ext_r i_GSM_r phi_GSMP_r w_GSM_r w_P_r


f_M_r = subs(d_x_M, [i_GSM   phi_GSMP   w_GSM   w_P   u_GSM   M_ext], ...
                      [i_GSM_r phi_GSMP_r w_GSM_r w_P_r u_GSM_r M_ext_r]);

% TODO                  
%i_GSM_r = solve(f_M_r(1), i_GSM_r);
%w_GMS_r = simplify(subs(solve(f_M_r(3), w_GSM_r)));
%w_P_r = solve(subs(f_M_r(4), w_P_r));


% Systemmatrizen des linearisierten Systems
A = [diff(f_M, i_GSM) diff(f_M, phi_GSMP) diff(f_M, w_GSM) diff(f_M, w_P)];
A = simplify(subs(A, w_P, w_P_r));
bu = [diff(f_M, u_GSM)];
bd = [diff(f_M, M_ext)];
ct = [0 0 0 1];

% Eigenwerte der Systemmatrix A berechnen
paralist_1 = [u_GSM_r M_ext_r L_GSM R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];

A_num = simplify(subs(A, paralist_1, paralist_2))

