clear all; close all; clc

%% Aufgabe 2.8.1 Abtastzeit bestimmen

V = 1.5;
xi = 0.9;
T = 0.5;
t_r = 1.2;
ue = 10;

Gs = tf([V], [T^2 2*xi*T 1]);

R_1 = tf([1], [1 0]);

omega_c = 1.5/t_r;
phi_soll = 70 - ue;

L_1 = Gs*R_1;

% Phasenreserve bei L_1(I*omega_c)
[re_L_1 im_L_1] = nyquist(L_1, omega_c);
phi_L_1 = atan (im_L_1/re_L_1) * 180/pi % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,
                                              % - 180 um in den richtigen
                                              % Quadranten zu kommen.
phi_dif = phi_soll - phi_L_1;

% Phase muss um phi_dif = 31.557 Grad angehoben werden
T_I = tan(phi_dif * pi /180)/omega_c;
R_2 = tf([T_I 1],1) * R_1;
L_2 = R_2*Gs;

% Phasenreserve bei L_2(I*omega_c)
[re_L_2 im_L_2] = nyquist(L_2, omega_c);
phi_L_2 = atan (im_L_2/re_L_2) * 180/pi

% Betrag korrigieren, mit dem Verstaerkungsfaktor
abs_L_2 = sqrt(re_L_2^2 + im_L_2^2) % V_R*abs(L_2(I*omega_c) = 1
V_R = 1/(abs_L_2)
R_3 = V_R*R_2;
L_3 = R_3*Gs;

%% D: Bodediagramm und ueberpruefen ob die Bedingungen erfuellt sind 

figure
line([omega_c omega_c], [25, -25])
hold on
bode(Gs, L_1, L_2, L_3)
%hold on
line([omega_c omega_c], [-90, -160])
line([400 600], [phi_soll-180,phi_soll-180])
grid on
title('2.1: Bode-Diagramm von G(s) und der offfene Regelkreise')
legend('G(s)', 'L1(s)','L2(s)','L3(s)', 'omega_c', 'phi soll');

L = L_3;
R = R_3;

T_ry = minreal(L/(1+L));

figure
step(T_ry)
% Ueberschwingung einzeichnen
line([0, 10], [1.10, 1.10], 'Color', 'r')
% tr so halbwegs einzeichnen, wie Abbildung 5.2.
a = 0.78; % Wendepunkt, vom Plot abgelesen (anklicken)
line([a-t_r/2, a+t_r/2], [0, 1], 'Color','k')
line([a-t_r/2, a-t_r/2], [0, 1], 'Color','g')
line([a+t_r/2, a+t_r/2], [0, 1], 'Color','g')
title('Sprungantwort des geschlossenen Kreises L3(s)')
legend('Try', 'ue', 'tr')
grid on

%%

t_rg = T*exp((acos(xi))/(tan(acos(xi))))

Ta_min = t_rg/10
Ta_max = t_rg/4

Ta = t_rg/40 % gewählt


%% Regler im q Bereich

Gz = c2d(Gs, Ta, 'zoh');
Gq = d2c(Gz, 'tustin');

Omega_c = 1.2/t_r;
phi_soll = 70 - ue;

Rq_1 = tf([1],[1,0]);
Lq_1 = Gq*Rq_1;

% Phasenreserve bei Ll_1(I*Omega_c)
[re_Ll_1 im_Ll_1] = nyquist(Lq_1, Omega_c);
phi_Lq_1 = atan (im_Ll_1/re_Ll_1) * 180/pi % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,
                                              % - 180 um in den richtigen
                                              % Quadranten zu kommen.

phi_dif = phi_soll - phi_Lq_1

% Phase um 24.3958 heben
Tq_I = tan(phi_dif * pi /180)/Omega_c;
Rq_2 = tf([Tq_I 1],1) * Rq_1;
Lq_2 = Rq_2*Gq;

% Phasenreserve bei L_2(I*omega_c)
[re_L_2 im_L_2] = nyquist(Lq_2, Omega_c);
phi_L1_2 = atan (im_L_2/re_L_2) * 180/pi


% Betrag korrigieren, mit dem Verstaerkungsfaktor
abs_Lq_2 = sqrt(re_L_2^2 + im_L_2^2) % V_R*abs(L_2(I*omega_c) = 1
Vq_R = 1/(abs_Lq_2)
Rq_3 = Vq_R*Rq_2;
Lq_3 = Rq_3*Gq;

Lq = Lq_3;
Rq = Rq_3;

%% Bode

figure
line([Omega_c Omega_c], [25, -25])
hold on
bode(Gq, Lq_1, Lq_2, Lq_3)
%hold on
line([Omega_c Omega_c], [-90+360, -160+360])
line([0.5 2], [phi_soll+180,phi_soll+180])
grid on
title('8.2: Bode-Diagramm von G#(q) und der offfene Regelkreise')
legend('G#(q)', 'L1(s)','L2(s)','L3(s)', 'omega_c', 'phi soll');

%% Steps

% Gz

Rz1 = c2d(R, Ta, 'zoh');
Rz2 = c2d(Rq, Ta, 'tustin');

Lz1 = Gz*Rz1;
Lz2 = Gz*Rz2;

Try1 = Lz1 / (1+ Lz1);
Try2 = Lz2 / (1+ Lz2);

figure
step(Try1, Try2)
% Ueberschwingung einzeichnen
line([0, 10], [1.10, 1.10], 'Color', 'r')
grid on
legend ('Try1', 'Try2')
title('Abtastzeit = 0.0317s')

stepinfo(Try1)
