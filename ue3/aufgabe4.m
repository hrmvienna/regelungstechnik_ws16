clear all; close all; clc

%% Aufgabe 3.4: Stoerverhalten im Vergleich

Ta = 10e-3
Ta2 = 50e-3

% Kompensationsregler mit Ta=10ms aus Aufgabe 2.5

R_komp_num = [0.154978970147333,-0.460423308714992,0.457094921736261,-0.151642123722946];
R_komp_den = [1,-2.86473429951691,2.73404280146561,-0.869308501948704];

R_komp = tf(R_komp_num, R_komp_den, Ta);

% Kompensationsregler mit Ta=50ms aus Aufgabe 2.5

R_komp2_num = [0.129052836691455,-0.350640525783673,0.338341661652813,-0.115965884301370];
R_komp2_den = [1,-2.40425531914894,1.89723856948846,-0.492983250339520];

R_komp2 = tf(R_komp2_num, R_komp2_den, Ta2);

% PI-Zustandsregler aus Aufgabe 3.3

R_pi_kI = 0.0065;
R_pi_kP = 0.1488;
R_pi_kx = [-0.6050   -0.7335    0.3827];

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

x_r = [i_GSM_r; phi_GSMP_r; w_GSM_r; w_P_r];
