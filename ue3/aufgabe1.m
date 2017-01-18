clear all; close all; clc

%% Aufgabe 3.1.1: Zustandsreglerentwurf

% Parameter
Ta = 10e-03; % 10ms

% Lin. red. Modell aus Aufgabe 2.3
A = double([0, 1, -1; -3411/62, -54163/28520, 1/1240; ...
    6822/325, 1/3250, 54117/74750 - (7*1130385^(1/2))/7475]);
bu = double([0; 12500/713; 0]);
bd = double([0; 0; -400/13]);
ct = [0 0 1];
d = 0;
sys = ss(A,[bu, bd], ct, d)

% System abtasten
dsys = c2d(sys, Ta)

phi = dsys.A;
gamma = dsys.B(:,1); % Spalte 1

% Erreichbarkeitsmatrix
R = ctrb(phi, gamma);

if rank(R) == 3
   sprintf('System ist vollstaendig erreichbar')
else
    sprintf('System ist nicht vollstaendig erreichbar')
end

%% Aufgabe 3.1.2: diskreter Zustandsregler

syms phi_GSMP w_GSM w_P uk rk g
kt = sym('kt', [1 3]);

xk = [phi_GSMP w_GSM w_P].';
uk = kt*xk + g*rk;

dxk = phi*xk + gamma*uk;
yk = ct*xk + d*uk;

% Ueberpruefung, dxk == dxk2, Skript Kapitel 8.1
phi_g = phi + gamma*kt;
dxk2 = phi_g*xk + gamma*g*rk;

% Eigenwerte der Dynamikmatrix
eigen_phi = eig(phi);

% Polvorgabe, -1/2 wurde gewaehlt, da es kleiner als die Eigenwerte von Phi
% ist, d.h. schnellere Dynamik der gewaehlten Pole als von Phi
lambda0 = -1/2;
P = [lambda0, lambda0, lambda0];
kt = -acker(phi, gamma, P) % Vorzeichen, Hinweis in Beispiel 8.1
kt1 = kt(1); kt2 = kt(2); kt3 = kt(3);

% Gleichung (8.39)
g = 1 / ((ct + d*kt)*inv(eye(3) - phi - gamma*kt)*gamma + d)

phi_g = double(subs(phi_g));
% Streckenuebertragungsfunktion
sys_g = ss(phi_g, gamma*g, (ct+ d*kt), d*g, Ta);
Gz = tf(sys_g)

%% Reglerentwurf
%% A: Anforderungen aus Aufgabe 2.5

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

omega_c = 1.2/t_r
phi_soll = 70 - u_e

% Da e_inf = 0 auf die Eingangsfolge (r^k) = (1^k), wird mind. ein
% Integrator benoetigt.

%% C: Reglerentwurf: zuerste Phase und Verstaerkung der bekannten Terme 

Gq = d2c(Gz, 'tustin');

% Regler aus den bekannten Termen
Rq_1 = tf([1], [1 0]);
Lq_1 = minreal(Rq_1*Gq);

% Phasenreserve bei Lz_1(I*omega_c)
[re_Lq_1 im_Lq_1] = nyquist(Lq_1, omega_c);
phi_Lq_1 = atan (im_Lq_1/re_Lq_1) * 180/pi; % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,

phi_dif = phi_soll - phi_Lq_1;

% Phase muss um 41.9121 gehoben werden
% mithilfe des Terms (1 + s*T_I)
% arctan(omega_c * T_I) = phi_dif * pi/180
T_I = tan(phi_dif * pi /180)/omega_c;
Rq_2 = tf([T_I 1], 1);
Lq_2 = minreal(Rq_2*Lq_1);

% Phasenreserve bei Lz_1(I*omega_c)
[re_Lq_2 im_Lq_2] = nyquist(Lq_2, omega_c);
phi_Lq_2 = atan (im_Lq_2/re_Lq_2) * 180/pi; % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,

% Betrag korregieren, mit dem Verstaerkungsfaktor
abs_Lq_2 = sqrt(re_Lq_2^2 + im_Lq_2^2); % V_R*abs(L_3(I*omega_c) = 1

V_R = 1/(abs_Lq_2);
Lq_3 = minreal(V_R*Lq_2);

% Gesamtregler
Rq = Rq_1*Rq_2*V_R;
Rz = c2d(Rq, Ta, 'tustin');

Lz = Rz*Gz;
Lz2 = c2d(Lq_3, Ta, 'tustin');

% residue von Rz_PI
[I,p1,P] = residue(cell2mat(Rz.Numerator), cell2mat(Rz.Denominator));


