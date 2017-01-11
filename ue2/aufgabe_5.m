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

%% Reglerentwurf
%% A: Streckenuebertragungsfunktion
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


%% B: Kenngroessen berechnen: Durchtrittsfrequenz und Phasenreserve

% Quadratischer Term der Strecke aus den konj.compl. Polstellen errechnen
pole_sys = pole(Gq);
p_con1 = pole_sys(1);
p_con2 = pole_sys(2);
syms x
coeff_p = coeffs(simplify(((x - p_con1))*(x - (p_con2))), x);
c1 = coeff_p(2);
c2 = coeff_p(1);
xi2T = c1/c2;
Tsq  = 1/c2;
T = double(sqrt(Tsq))
xi = double(c1/(2*c2*T))

% Kenngroessen berechnen
omega_c = 1.2/t_r
phi_soll = 70 - u_e

% Da e_inf = 0 auf die Eingangsfolge (r^k) = (1^k), wird mind. ein
% Integrator benoetigt.

%% C: Reglerentwurf: zuerste Phase und Verstaerkung der bekannten Terme 

% Regler aus den bekannten Termen
Rq_kompterm = tf([(T^2) (2*xi*T) 1], [1]);
V_1 = double(c2);
Rq_1 = V_1*tf([1], [1 0]);
Lq_1 = minreal(Rq_1*Rq_kompterm*Gq);

% % Realisierungpole waehlen
T_Real = 7;
Rq_real = tf([1], [1 T_Real]);
Rq_2 = Rq_real*Rq_real;
Lq_2 = Rq_2*Lq_1;

% Phasenreserve bei Lq_2(I*omega_c)
[re_Lq_2 im_Lq_2] = nyquist(Lq_2, omega_c);
phi_Lq_2 = atan (im_Lq_2/re_Lq_2) * 180/pi; % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,

phi_dif = phi_soll - phi_Lq_2

% Phase muss um 41.9121 gehoben werden
% mithilfe des Terms (1 + s*T_I)
% arctan(omega_c * T_I) = phi_dif * pi/180
T_I = tan(phi_dif * pi /180)/omega_c;
Rq_3 = tf([T_I 1], 1);
Lq_3 = Rq_3*Lq_2;

% Phasenreserve bei Lq_3(I*omega_c)
[re_Lq_3 im_Lq_3] = nyquist(Lq_3, omega_c);
phi_Lq_3 = atan (im_Lq_3/re_Lq_3) * 180/pi;

% Betrag korregieren, mit dem Verstaerkungsfaktor
abs_Lq_3 = sqrt(re_Lq_3^2 + im_Lq_3^2) % V_R*abs(L_3(I*omega_c) = 1
V_R = 1/(abs_Lq_3)
Rq_4 = V_R;
Lq_4 = Rq_4*Lq_3;

% Regler im q-Bereich in PI- und Kompensationsteil teilen
Rq_komp = Rq_kompterm*Rq_2
Rq_PI = Rq_1*Rq_3*Rq_4

% Gesamtregler
Rq = Rq_PI*Rq_komp;

%% D: Bodediagramm und ueberpruefen ob die Bedingungen erfuellt sind 

figure
line([omega_c omega_c], [5, -5], 'Color','k');
hold on
h = bodeplot(Gq, Lq_1, Lq_2, Lq_3, Lq_4);
p = getoptions(h); % return plotoptions for bode plot
p.PhaseMatching = 'on';
p.PhaseMatchingValue = 0;
setoptions(h,p);
%bode(Gq, Lq_1, Lq_2, Lq_3)
line([omega_c omega_c], [-100, -120], 'Color','k');
grid on
title('2.1: Bode-Diagramm von G#(q) und der offfene Regelkreise')
legend('G#(q)', 'L1(q)','L2(q)','L3(q)','L4(q)', 'omega_c');

%% E: Sprungantwort des geschlossenen Kreises

% Geschlossener Kreis
T_ry = Lq_4 / (1 + Lq_4);

figure
step(T_ry)
% Ueberschwingung einzeichnen
line([0, 7], [1, 1], 'Color', 'r')
% tr so halbwegs einzeichnen, wie Abbildung 5.2.
a = 0.75; % Wendepunkt, vom Plot abgelesen (anklicken)
line([a-t_r/2, a+t_r/2], [0, 1], 'Color','k')
line([a-t_r/2, a-t_r/2], [0, 1], 'Color','g')
line([a+t_r/2, a+t_r/2], [0, 1], 'Color','g')
title('Sprungantwort des geschlossenen Kreises L4(s)')
legend('Try', 'ue', 'tr')
grid on

%% F: Stellegroessenanforderungen

% Stellgroesse darf nur im Bereich von 0 bis 12 V liegen.
% Ueberpruefung in der Simulation aufgabe_6_sim.slx

%% Regler im z-Bereich

Rz = c2d(Rq, Ta, 'tustin');
Rz_komp = c2d(Rq_komp, Ta, 'tustin');
Rz_PI = c2d(Rq_PI, Ta, 'tustin');
% Ueberpruefung: Rz == Rz_PI*Rz_komp => Ja

% residue von Rz_PI
[I,p1,P] = residue(cell2mat(Rz_PI.Numerator), cell2mat(Rz_PI.Denominator));

% Proportional und Integralglied
Rz_p = tf([P],[1], Ta);
Rz_i = tf([I],[1 -p1], Ta);
Rz_p_i = Rz_p + Rz_i;

% Ueberpruefung: Rz_p_i == Rz_PI => Ja

% Gesamtregler im Z Bereich
Rz_gesamt = Rz_p_i*Rz_komp;

