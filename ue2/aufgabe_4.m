clear all; close all; clc

%% Aufgabe 2.4.1

syms u_GSM M_ext i_GSM phi_GSMP w_GSM w_P M_GSM M_kopp

u = [u_GSM M_ext].'
y1 = [i_GSM phi_GSMP w_GSM w_P].'
y2 = [M_GSM M_kopp].'

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

% TODO Angaben anpassen, falls sich was in Aufgabe 2.2/2.3 aendert
% Ruhelagen
i_GSM_r = (u_GSM_r + (k_GSM*(R_GSM*d_vP + R_GSM*d_vGSM - ... 
    (k_GSM^4 + R_GSM^2*d_vP^2 + R_GSM^2*d_vGSM^2 - ...
    4*M_ext_r*R_GSM^2*d_qP - 4*R_GSM^2*d_cP*d_qP - ...
    4*R_GSM^2*d_cGSM*d_qP + 2*R_GSM^2*d_vP*d_vGSM + ...
    2*R_GSM*d_vP*k_GSM^2 + 2*R_GSM*d_vGSM*k_GSM^2 + ...
    4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + k_GSM^2))/(2*R_GSM*d_qP))/R_GSM;

phi_GSMP_r = ((k_GSM^2 + R_GSM*d_vGSM)*(R_GSM*d_vP + R_GSM*d_vGSM - ...
    (k_GSM^4 + R_GSM^2*d_vP^2 + R_GSM^2*d_vGSM^2 - ...
    4*M_ext_r*R_GSM^2*d_qP - 4*R_GSM^2*d_cP*d_qP - ...
    4*R_GSM^2*d_cGSM*d_qP + 2*R_GSM^2*d_vP*d_vGSM + ...
    2*R_GSM*d_vP*k_GSM^2 + 2*R_GSM*d_vGSM*k_GSM^2 + ...
    4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + k_GSM^2))/(2*R_GSM^2*c_GSMP*d_qP)...
    - (R_GSM*d_cGSM - k_GSM*u_GSM_r)/(R_GSM*c_GSMP);

w_GSM_r = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

w_P_r = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

w_P_r_red = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

% Anfangsbedingungen fuers Nicht-Lineare System
i_GSM_0 = i_GSM_r;
phi_GSMP_0 = phi_GSMP_r;
w_GSM_0 = w_GSM_r;
w_P_0 = w_P_r;

% Linearisiertes Modell
A_lin = [ -R_GSM/L_GSM, 0, -k_GSM/L_GSM, 0; ...
    0, 0, 1, -1; ...
    k_GSM/J_GSM, -c_GSMP/J_GSM, -(d_GSMP + d_vGSM)/J_GSM, d_GSMP/J_GSM; ...
    0, c_GSMP/J_P, d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r)/J_P];

bu_lin = [1/L_GSM; 0; 0; 0];
bd_lin = [0; 0; 0; -1/J_P];
b_lin = [bu_lin bd_lin];
c_lin = [0 0 0 1];
d_lin = [0 0];

% Reduziertes Modell
A_red = [0,1,-1; ...
-c_GSMP/J_GSM, -(k_GSM^2/R_GSM + d_GSMP + d_vGSM)/J_GSM, d_GSMP/J_GSM; ...
c_GSMP/J_P, d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r_red)/J_P];

bu_red = [0; k_GSM/(J_GSM*R_GSM); 0];
bd_red = [0; 0; -1/J_P];
b_red = [bu_red bd_red];
c_red = [0 0 1];
d_red = d_lin;

