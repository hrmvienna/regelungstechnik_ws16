clear all; close all; clc

%% Antwort zu 2.6
% Alle drei Systeme weisen kein Überschwingen auf, die Anstiegszeit wird 
% eingehalten. Es bleibt keine Regelabweichung und die Stellgröße bleibt 
% im Bereich [0, 12] V (Siehe aufgabe_6.fig).
% Der Ausgang des nicht-Linearisiertem System steigt etwas spaeter an und 
% schwingt noch leicht nach aber nicht über. 
% Die beiden linearisierten Systeme schwingen minimal vor dem Erreichen des
% Endwertes.

%% Parameter und Ruhelage

% Zustandsgroessen
syms delta_phi_GSMP delta_w_GSM delta_w_P w_P_r u_GSM_r M_ext_r
% GSM
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
% P
syms J_P d_cP d_vP d_qP
% GSMP
syms c_GSMP d_GSMP phi_GSMP

syms u_GSM M_ext i_GSM phi_GSMP w_GSM w_P M_GSM M_kopp

u = [u_GSM M_ext].';
y1 = [i_GSM phi_GSMP w_GSM w_P].';
y2 = [M_GSM M_kopp].';

Ta = 50e-03;

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

%% Reduzierte System
A_red = [             0,                                        1,                               -1;...
 -c_GSMP/J_GSM, -(k_GSM^2/R_GSM + d_GSMP + d_vGSM)/J_GSM,                        d_GSMP/J_GSM;...
    c_GSMP/J_P,                               d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r)/J_P];
bu_red = [0; k_GSM/(J_GSM*R_GSM); 0];
bd_red = [0; 0; -1/J_P];
ct_red = [0 0 1];
d_red = 0;

sys_red = ss(A_red, [bu_red, bd_red], ct_red, d_red);
tf_sys_red = tf(sys_red);

% Uebertragungsfunktionen T_ry und T_dy, wobei r der Eingang und d die
% Stoerung ist
Gu_red = tf_sys_red(1);
Gd_red = tf_sys_red(2);

% Zeitdiskrete Uebertragungsfunktion T_ry und T_dy
Gz_red = c2d(Gu_red, Ta, 'zoh');
Gzd_red = c2d(Gd_red, Ta, 'zoh');
Gq_red = d2c(Gz_red, 'tustin');
Gqd_red = d2c(Gzd_red, 'tustin');

%% Linearisiertes Modell
A_lin = [ -R_GSM/L_GSM, 0, -k_GSM/L_GSM, 0; ...
    0, 0, 1, -1; ...
    k_GSM/J_GSM, -c_GSMP/J_GSM, -(d_GSMP + d_vGSM)/J_GSM, d_GSMP/J_GSM; ...
    0, c_GSMP/J_P, d_GSMP/J_P, -(d_GSMP + d_vP + 2*d_qP*w_P_r)/J_P];

bu_lin = [1/L_GSM; 0; 0; 0];
bd_lin = [0; 0; 0; -1/J_P];
b_lin = [bu_lin bd_lin];
c_lin = [0 0 0 1];
d_lin = [0 0];

% Uebertragungsfunktion lin. System
Gs_lin = ss(A_lin, b_lin, c_lin, d_lin);
tf_sys_lin = tf(Gs_lin);

% Uebertragungsfunktionen T_ry und T_dy, wobei r der Eingang und d die
% Stoerung ist
Gu_lin = tf_sys_lin(1);
Gd_lin = tf_sys_lin(2);

% Zeitdiskrete Uebertragungsfunktion T_ry und T_dy
Gz_lin = c2d(Gu_lin, Ta, 'zoh');
Gzd_lin = c2d(Gd_lin, Ta, 'zoh');
Gq_lin = d2c(Gz_lin, 'tustin');
Gqd_lin = d2c(Gzd_lin, 'tustin');

%% PI-Komp Regler aus Aufgabe 2.5
syms z

Ta = 50e-03;
Rz_num = [0.129052836691455 -0.350640525783673 0.338341661652813 -0.115965884301370];
Rz_den = [1 -2.40425531914894 1.89723856948846 -0.492983250339520];
Rz = tf(Rz_num, Rz_den, Ta);
Rq = d2c(Rz, 'tustin');

% Regler mit t_r = 4
num = [0.0114235121008542 -0.0302399436766243 0.0285518173838854 -0.00952261914285287];
den = [1 -2.40425531914894 1.89723856948846 -0.492983250339520];
Rz_tr4 = tf(num, den, Ta);

% Regler mit t_r = 0.5, realisierungspol = 15
num = [0.856133395773891 -2.32771511276413 2.24730725704585 -0.770779483086807];
den = [1 -1.90909090909091 1.11570247933884 -0.206611570247934];
Rz_tr05 = tf(num, den, Ta);



% Zusaetzliche Parameter
te = 1; % Start Eingangssprung
ts = 0.5; % Start Stoersprung
tr = 5; % Start Rampe
Usinus = 1; % sin amp
we = 2; % sin omega

% figure
% step((Rq*Gq)/(1 + Rq*Gq))
% grid on
