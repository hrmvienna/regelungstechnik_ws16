clear all; close all; clc

%% Aufgabe 2.4.1

syms u_GSM M_ext i_GSM phi_GSMP w_GSM w_P M_GSM M_kopp

u = [u_GSM M_ext].'
y1 = [i_GSM phi_GSMP w_GSM w_P].'
y2 = [M_GSM M_kopp].'

u_GSM_r = 5.6;
M_ext_r = 0;
i_GSM_0 = 0;
phi_GSMP_0 = 0;
w_GSM_0 = 0;
w_P_0 = 0;

% Parameter
L_GSM   = 1.4;
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
