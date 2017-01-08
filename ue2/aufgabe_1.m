clear all   % Es werden alle Variablen geloescht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurueckgesetzt

%% Aufgabe 2.1.1 PI Regler

%
% Vorgangsweise beim Reglerentwurf nach dem Frequenzkennlinienverfahren
% ---------------------------------------------------------------------
%
% (A) Zu einer gegebenen Streckenuebertragungsfunktion G(s) muessen die Kenngroessen des
% Einschwingverhaltens des geschlossenen Kreises (t r , M oder ue und e_inf ) spezifiziert
% werden.
%
% (B) Die Kenngroessen t r , M oder ue und e_inf werden mithilfe der Beziehungen (5.2), (5.3)
% und (5.9) in Vorgaben an den Frequenzgang des offenen Kreises L(I omega) uebersetzt.
%
% (C) Ein Regler R(s) muss so gewaehlt werden, dass der geschlossene Kreis BIBO-stabil
% ist und die Forderungen von (B) erfuellt werden. Erfuellt die Uebertragungsfunktion
% des offenen Kreises L(s) = R(s)G(s) die Bedingungen von Satz 4.6, dann kann die
% Stabilitaet des geschlossenen Kreises anhand der Phasenreserve Omega beurteilt werden,
% anderenfalls muss man das Nyquistkriteriums von Satz 4.5 anwenden.
%
% (D) Um ein kriechendes Einlaufen der Sprungantwort in den stationaeren Endwert zu
% vermeiden, soll in (C) der Regler R(s) so entworfen werden, dass ca. 1 Dekade
% um die Durchtrittsfrequenz Omega_C die Betragskennlinie von L(s) mit mindestens 20
% dB/Dekade abfaellt.
%
% (E) Die Qualitaet des Entwurfes ist immer durch Simulation zu ueberpruefen, insbesondere
% auch deshalb, weil das Verfahren sich auf empirische Formeln stuetzt. Sind die
% Ergebnisse nicht zufriedenstellend, dann muss man sich die Frage stellen, ob die
% Anforderungen von (A) ueberhaupt prinzipiell erfuellbar sind, oder ob ein anderer
% Regler R(s) von (C) die Situation verbessern wuerde.
%
% (F) Die Begrenzung der Stellgroesse u(t), die bei jedem technisch relevanten Prozess vor-
% handen ist, kann im Rahmen dieses einfachen Entwurfsverfahrens nicht systematisch
% beruecksichtigt werden. Sollte sich bei der Simulation herausstellen, dass man zu viel
% Stellgroesse benoetigt, dann muss man die Anforderungen in (A) entsprechend den
% Ueberlegungen von Abschnitt 4.3.1 veraendern, also die Anstiegszeit t r vergroessern. Im
% Rahmen einer Fuehrungsregelung sollte auf keinen Fall ein Sprung sondern immer
% ein hinreichend glattes Signal als Fuehrungsgroesse verwendet werden (man wiederhole
% dazu auch die ueberlegungen von Abschnitt 4.3.2).

%% A: Streckenuebertragungsfunktion

syms s rho r_d T_t V xi T

cc = 9.959e05; % = T^2, herausgehobener Term, zwecks Umformen

% Uebertragungsfunktion

% Zaehlerpolynom
z_L = V;
% Zaehlerpolynom Stoerfunktion
z_Ld = (857.1/cc) * s + (1.616e06/cc);
% Nennerpolynom
n_L = 1 + 2*xi*(s*T) + (s*T)^2;
% Uebertragungsfunktion
G = z_L / n_L

% Stoerfunktion
G_d = z_Ld / n_L

% Anforderungen
t_r=3e-3
u_e=5
e_a_inf = 0

%% B: Kenngroessen berechnen: Durchtrittsfrequenz und Phasenreserve
% Der Regler muss mind. eine Polstelle bei s = 0 haben, da der Regelfehler
% fuer die Sprungantwort e_inf = 0. (5.10, Seite 125)

% phi - phasenreserve
% u_e - Prozentuelles Ueberschwingen

% (5.2) omega_c * tr ~~ 1.5
omega_c = 1.5/t_r

% (5.3) phi[grad] + u_e[%] = 70
phi_soll = 70 - u_e

% (5.9) e_inf = lim[s->0]((s*()s^p*n_L(s)) / (s^p*n_L(s) +
% V*z_L(s)*e^(-sT_t)))*r_d(s)

e_inf = limit(s * ((s^rho*n_L) / (s^rho*n_L + V*z_L*exp(-s*T_t))) * r_d, s, 0);

% Reglerabweichung fuer die Sprungantwort, _s = sprung
r_d1 = 1/s;
% rho = 0
e_inf1_s = subs(e_inf, [r_d, rho], [r_d1, 0])
% rho = 1
e_inf2_s = subs(e_inf, [r_d, rho], [r_d1, 1])

% Reglerabweichung fuer die Rampenantwort, -r = rampe
r_d2 = 1/(s^2);
% rho = 0
e_inf1_r = subs(e_inf, [r_d, rho], [r_d2, 0])
% rho = 1
e_inf2_r = subs(e_inf, [r_d, rho], [r_d2, 1])

% Schlussfolgerung:
% e_inf2_s => V != 0 und
% e_inf2_r => 1/V^2 = 1e-3 fuer Aufgabe 2.1.3
% !!! e_inf2_r sollte aber 1/V sein, nicht 1/V^2

%% C Reglerentwurf: zuerste Phase und Verstaerkung der bekannten Terme 
% berechnen und danach diese Korrekieren.
% R_1 = 1/s
T = 1.002e-3;
xi = 0.802;
V = 1.377;

% G(s) = V/(T^2*s^2 + 2*xi*T*s + 1)
V_G = V;

% Uebertragungsfunktion der Strecke, G(s)[Gs] und Gd(s)
b = [V_G];
a = [(T^2) 2*xi*T 1];
c = [857.1/cc 1.616e06/cc];
Gs = tf(b, a)
Gd = tf(c, a)
Gsys = [Gs, Gd]; % Uebertragungsfunktioen Eingang/Stoerung zum Ausgang

% Regler mit allen bisher bekannten Termen, R_1(s) = 1/s
R_1 = tf(1,[1 0])

% Uebertragungsfunktion des offenen Regelkreises L_1(s) = R_1(s)*G(s)
% wobei fuer R(s) ein PI Regler verwendet wurde. (Siehe [5.1 - PI-Reglerentwurf]
L_1 = R_1*Gs

% Phasenreserve bei L_1(I*omega_c)
[re_L_1 im_L_1] = nyquist(L_1, omega_c);
phi_L_1 = atan (im_L_1/re_L_1) * 180/pi % phi = arctan(Im/re)[rad], [degree] = [rad]*180/pi,
                                              % - 180 um in den richtigen
                                              % Quadranten zu kommen.
phi_dif = phi_soll - phi_L_1

% Phase muss um phi_dif = 22.0143 Grad angehoben werden
% mithilfe des Terms (1 + s*T_2)
% arctan(omega_c * T_2) = phi_L_1 * pi/180
T_2 = tan(phi_dif * pi /180)/omega_c;
R_2 = tf([T_2 1],1) * R_1;
L_2 = R_2*Gs;

% Phasenreserve bei L_2(I*omega_c)
[re_L_2 im_L_2] = nyquist(L_2, omega_c);
phi_L_2 = atan (im_L_2/re_L_2) * 180/pi


% Betrag korregieren, mit dem Verstaerkungsfaktor
abs_L_2 = sqrt(re_L_2^2 + im_L_2^2) % V_R*abs(L_2(I*omega_c) = 1
V_R = 1/(abs_L_2)
R_3 = V_R*R_2;
L_3 = R_3*Gs;

%% D: Bodediagramm und ueberpruefen ob die Bedingungen erfuellt sind 

figure
line([omega_c omega_c], [25, -150])
hold on
bode(Gs, L_1, L_2, L_3)
%hold on
line([omega_c omega_c], [-90, -160])
line([400 600], [phi_soll-180,phi_soll-180])
grid on
title('2.1: Bode-Diagramm von G(s) und der offfene Regelkreise')
legend('G(s)', 'L1(s)','L2(s)','L3(s)', 'omega_c', 'phi soll');

%% D: Sprungantwort des geschlossenen Kreises

% Geschlossener Kreis
T_ry = L_3 / (1 + L_3);

figure
step(T_ry)
% Ueberschwingung einzeichnen
line([0, 1], [1.05, 1.05], 'Color', 'r')
% tr so halbwegs einzeichnen, wie Abbildung 5.2.
a = 0.002; % Wendepunkt, vom Plot abgelesen (anklicken)
line([a-t_r/2, a+t_r/2], [0, 1], 'Color','k')
line([a-t_r/2, a-t_r/2], [0, 1], 'Color','g')
line([a+t_r/2, a+t_r/2], [0, 1], 'Color','g')
title('Sprungantwort des geschlossenen Kreises L3(s)')
legend('Try', 'ue', 'tr')
grid on

%% F: Stellegroessenanforderungen

% keine

%% Aufgabe 2.1.2 PI Regler

% Ist das FKL-Verfahren ein exaktes Entwurfsverfahren, wenn die Strecke ein
% Verzoegerungsglied 2-ter Ordnung ist?

% Die Strecke ist ein Modell/Abbild der Wirklichkeit und ist dahe nicht
% 100%ig exakt. Des Weiteren baut das FKL-Verfahren auf numerisch
% errechneten Werten (Verhaeltnis omega_c zu t_r, Ueberschwingen) auf, und
% enthaelt somit Ungenauigkeiten, somit ist das Verfahren nicht exakt.

% TODO?

%% Aufgabe 2.1.3 PID Regler - Aehnlich: Lead-Lag-Reglerentwurf

% A: Streckenuebertragungsfunktion
% gleich wie Aufgabe 2.1
% B: Kenngroessen berechnen
% t_r und ue gleich wie Aufgabe 2.1

% Bleibende Regelabweichung fuer die Rampenantwort
e_inf_rampe = 1e-3;
V_L = 1/e_inf_rampe;

% V_R_pid*V_G = V_L;
V_R_pid = V_L/V_G;

%% C: Reglerentwurf: PID Regler: R(s) = Vp * (1 + TI*s)(1 + TD*s)/s(1 + TR*s)

% Zeitkonstante des Intagralterms
T_I = 1e-3;

Rr_1 = V_R_pid*tf([T_I 1],[1 0]);
Ll_1 = Rr_1*Gs

% Phasenreserve bei Ll_1(I*omega_c)
[re_Ll_1 im_Ll_1] = nyquist(Ll_1, omega_c);
phi_Ll_1 = atan (im_Ll_1/re_Ll_1) * 180/pi
abs_L1_1 = sqrt(re_Ll_1^2 + im_Ll_1^2)

% => Betrag und Phase muessen gesenkt werden
% => Lead-Glied nicht erforderlich, da die Phase schon hoch genug ist

% Lag Glied entwerfen: (Seite 131, 5.31)
delta_a = 1/abs_L1_1;
delta_phi = (phi_soll - phi_Ll_1) * pi/180; % [rad]

T_lag = (delta_a*sqrt(1 + tan(delta_phi)^2) - 1)/(omega_c*tan(delta_phi));
eta_lag = (omega_c*T_lag - tan(delta_phi))/(omega_c*T_lag* (1 + omega_c*T_lag*tan(delta_phi)));

% Lag Regler - senkt Betrag und Phase
R_lag = tf([T_lag 1],[(T_lag*eta_lag) 1])

% offener Kreis
Ll_2 = Rr_1*R_lag*Gs

% Phasenreserve bei Ll_2(I*omega_c)
[re_Ll_2 im_Ll_2] = nyquist(Ll_2, omega_c);
phi_Ll_2 = atan (im_Ll_2/re_Ll_2) * 180/pi
abs_L1_2 = sqrt(re_Ll_2^2 + im_Ll_2^2)

figure
line([omega_c omega_c], [25, -50])
hold on
bode(Ll_1, Ll_2)
%hold on
line([omega_c omega_c], [-90, -160])
line([400 600], [phi_soll-180,phi_soll-180])
title('Bode-Diagramm fuer den offenen Kreis L1(s) und L2(s) mit PID-Regler')
legend('L1(s)', 'L2(s)')
grid on

%% Aufgabe 2.1.4: Geschlossener Kreis des lin. Systems als Blockschaltbild

% Blockschaltbild
% siehe aufgabe_1_4.slx

% Fuehrungsuebertragungsfunktion Try
T_ry_pi = L_3 / (1 + L_3)
T_ry_pid = Ll_2 / (1 + Ll_2)

% Stoeruebertragungsfunktion Tdy
T_dy_pi = Gd / (1 + L_3)
T_dy_pid = Gd / (1 + Ll_2)

% Sprungantworten
figure
step(T_ry_pi, T_dy_pi, T_ry_pid, T_dy_pid)
title('Sprungantworten fuer die Fuehrungs/Stoeruebertragungsfunktion')
legend('Try PI', 'Tdy PI', 'Try PID', 'Tdy PID')
grid on

%% Aufgabe 2.1.5: Simulieren der implementierten Regler

% Parameter aus Esys_Parameter.m von Aufgabe 1.4:
R1      = 625;    % Widerstand in Ohm
R2      = 3500;   % Widerstand in Ohm
K       = 3;      % Spannungsverstaerkung
C1_ref  = 1e-6;   % Referenz-Kapazitaetswert auf Kennlinie von C1 in F
UC1_ref = -10;    % Referenz-Spannungswert auf Kennlinie von C1 in V
kC1     = 800e-9; % Steigung der Kennlinie von C1 in F/V
C2      = 1e-6;   % Kapazitaet in F
UeR = 5;
UsR = 5;
te = 1;      % Einschaltzeitpunkt in Sekunden
Ue = 7;      % Endwert des Sprunges in Volt
Usinus = 0.1; % Amplitude des Sinus in Volt
we     = 10;  % Winkelfrequenz des Sinus in rad/s
ts = 3;     % Einschaltzeitpunkt in Sekunden
Us = 1;     % Endwert des Sprunges in Volt
UC1R = -(K*R2*UeR - R1*UsR - 2*R2*UeR + K*R1*UsR + K*R2*UsR)/(R1 + 2*R2);
UC2R = (R2*UeR + R1*UsR + R2*UsR)/(R1 + 2*R2);
UaR  = K*UC2R;
dQ1  = C1_ref + UC1R*kC1 - UC1_ref*kC1; % Erste Ableitung von Q1 nach UC1
ddQ1 = kC1; % Zweite Ableitung von Q1 nach UC1
C1 = dQ1;

A11 = -(R1 + R2)/(C1*R1*R2);
A12 = -(K*R1 - R1 + K*R2)/(C1*R1*R2);
A21 = 1/(C2*R2);
A22 = (K - 2)/(C2*R2);

A  = [A11,A12;A21,A22];
bu = [1/(C1*R1);0];
bd = [0;1/(C2*R2)];
c  = [0, K];
du = [0];
dd = [0];

%% Aufgabe 2.1.6: Messrauschen

% TODO: Verhalten mit und ohne Messrauschen beschreiben

%% Aufgabe 2.1.7: Steigungsbegraenzung

% TODO: Verhalten mit Steigungsbegrenzung untersuchen.


