clear all; close all; clc

%% Aufgabe 2.5: Kompensationsreglerentwurf im q-Bereich des linearisierten
%  Systems im Arbeitspunkt aus Aufgabe 2.3

% Zustandsgroessen
syms delta_phi_GSMP delta_w_GSM delta_w_P w_P_r u_GSM_r M_ext_r
% GSM
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
% P
syms J_P d_cP d_vP d_qP
% GSMP
syms c_GSMP d_GSMP phi_GSMP

% Parameterliste
paralist_1 = [u_GSM_r M_ext_r L_GSM R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];

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
Gu = tf_sys(1)
Gd = tf_sys(2)

%% Reglerentwurf
% Sollsprung delta_r = 20 rad s^-1
delta_r = 20;

% Bleibende Reglerabweichung fuer Eingangssprung r^k = (1^k) 
e_inf_s = 0;
% Antiegszeit
t_r = 1;
% Ueberschwingen
u_e = 0;
% Stellgroessenbeschraenkung
u_gsm_min = 0;
u_gsm_max = 12;

% Quadratischer Term des Reglers aus den konj.compl. Polstellen errechnen
pole_sys = pole(sys);
p_con1 = pole_sys(1);
p_con2 = pole_sys(2);
syms x
poly = simplify(((x - p_con1))*(x - (p_con2)));
xi2T = 1626954853497643/1125899906842624;
Tsq  = 382699032744171196714245354617545/5070602400912917605986812821504;
T = sqrt(Tsq)
xi = xi2T /(2*T)

% Regler aus den bekannten Termen
Rq_1 = tf([(T^2) (2*xi*T) 1], [1 0])
