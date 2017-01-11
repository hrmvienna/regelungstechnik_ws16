clear all; close all; clc

%% asd

%  Systems im Arbeitspunkt aus Aufgabe 2.3

% Zustandsgroessen
syms delta_phi_GSMP delta_w_GSM delta_w_P w_P_r u_GSM_r M_ext_r
% GSM
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
% P
syms J_P d_cP d_vP d_qP
% GSMP
syms c_GSMP d_GSMP phi_GSMP
Ta = 50e-03;

% Parameterliste
paralist_1 = [u_GSM_r M_ext_r L_GSM    R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4e-3   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];

% TODO: evtl anpassen, falls sich Aufgabe 2.3 aendert
% Linearisiertes System aus Aufgabe 2.3
w_P_r = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - 4*R_GSM^2*d_cP*d_qP - ...
    4*R_GSM^2*d_cGSM*d_qP + 2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

A = [             0,                                        1,                               -1;...
 -c_GSMP/J_GSM, -(k_GSM^2/R_GSM + d_GSMP + d_vGSM)/J_GSM,                        d_GSMP/J_GSM;...
    c_GSMP/J_P,                               d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r)/J_P];
bu = [0; k_GSM/(J_GSM*R_GSM); 0];
bd = [0; 0; -1/J_P];
ct = [0 0 1];
d = 0;

A = double(simplify(subs(A, paralist_1, paralist_2)));
bu = double(simplify(subs(bu, paralist_1, paralist_2)));
bd = double(simplify(subs(bd, paralist_1, paralist_2)));

sys = ss(A, [bu, bd], ct, d);
tf_sys = tf(sys);

% Uebertragungsfunktionen T_ry und T_dy, wobei r der Eingang und d die
% Stoerung ist
Gu = tf_sys(1);
Gd = tf_sys(2);

% Zeitdiskrete Uebertragungsfunktion T_ry und T_dy
Gz = c2d(Gu, Ta, 'zoh');
Gzd = c2d(Gd, Ta, 'zoh');
Gq = d2c(Gz, 'tustin');
Gqd = d2c(Gzd, 'tustin');

%% Systemmodelle aus Aufgabe 2.4.1

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

%% Aufgabe 2.6.1
syms z

Ta = 50e-03;
Rz_num = [0.129052836691455 -0.350640525783673 0.338341661652813 -0.115965884301370];
Rz_den = [1 -2.40425531914894 1.89723856948846 -0.492983250339520];
Rz = tf(Rz_num, Rz_den, Ta);
Rq = d2c(Rz, 'tustin');

Try_num = [0.000314580928418420 -0.00845210206034595 -1.68449864678759 37.1533505213221 838.487970229087 4130.14296061898 4649.67519742700 1508.98401881251 0];
Try_den = [1.00031458092842 29.4429060654243 333.480140175364 1850.59773547408 5385.57475499688 8337.36302265092 5914.06612575424 1508.98401881251 0];

Try = tf(Try_num, Try_den);

% Zusaetzliche Parameter
te = 1; % Start Eingangssprung
ts = 5; % Start Stoersprung
tr = 5; % Start Rampe
Usinus = 1; % sin amp
we = 2; % sin omega

figure
step((Rq*Gq)/(1 + Rq*Gq))
grid on