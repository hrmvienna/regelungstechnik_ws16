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

syms s rho r_d T_t V

% A:

% Uebertragungsfunktion

% Zaehlerpolynom
z_L = 1.371e06
% Nennerpolynom
n_L = s^2 + 1600 * s + 9.959e05
% Uebertragungsfunktion
G = z_L / n_L

% Stoerfunktion
G_d = (857.1 * s + 1.616e06) / (s^2 + 1600 * s + 9.959e05)

% Anforderungen
tr=3e-3
ue=5
e_a_inf = 0

% B:

% Errechne die Durchtrittsfrequenz, Phasenreserve

% (5.2) omega_c * tr ~~ 1.5
omega_c = 1.5/tr

% (5.3) phi[grad] + u_e[%] = 70
phi = 70 - ue

% (5.9) e_inf = lim[s->0]((s*()s^p*n_L(s)) / (s^p*n_L(s) +
% V*z_L(s)*e^(-sT_t)))*r_d(s)

% Sprungfunktion
r_d1 = 1/s

e_inf = limit(s * ((s^rho*n_L) / (s^rho*n_L + V*z_L*exp(-s*T_t))) * r_d, s, 0);

% rho = 0
e_inf1 = subs(e_inf, [r_d, rho], [r_d1, 0])
% rho = 1
e_inf2 = subs(e_inf, [r_d, rho], [r_d1, 1])


%%

% phi - phasenreserve
% u_e - Prozentuelles Ueberschwingen

% C:

% Uebertragungsfunktion der Strecke
b = [0 0 1.371e06]
a = [1 1600 9.959e05]
G_s_tf = tf(b, a)

% Uebertragungsfunktion der Regelstrecke L_1(s) = G(s)*R(s)
% wobei fuer R(s) ein PI Regler verwendet wurde. (Siehe [5.1 - PI-Reglerentwurf]
b2 = [0 0 0 1.371e06]
a2 = [1 1600 9.959e05 0]
L_1_s = tf(b2, a2)


bode(G_s_tf, L_1_s)
legend('G(s)', 'L1(s)');


% D:
% E:
% F:

%%

% Anforderungen an den offenen und geschlossenen Kreis
tr=3e-3;
ue=5;


%% Aufgabe 2.1.2 PI Regler

%Gesamtregler und offener Kreis
% R_PI = 
% L_PI = 
% bode(L_PI,'r');

%% Aufgabe 2.1.3 PID Regler

% Bleibende Regelabweichung
e_inf = 1e-3;

% Zeitkonstante des Intagralterms
TI = 1e-3;

% Gesamtregler und offener Kreis
% R_PID = 
% L_PID = 
% bode(L_PID,'g');

%% Aufgabe 2.1.4
% Kontrolle des Verhaltens des geschlossenen Kreise

T = 0.05;   % Simulationsdauer
figure;
subplot(3,1,1);hold on;grid on;
title('Sprungantwort Fuehrungsuebertragungsfunktion')

subplot(3,1,2);hold on;grid on;
title('Sprungantwort Stoeruebertragungsfunktion')

subplot(3,1,3);hold on;grid on;
title('Rampenantwort Stoeruebertragungsfunktion')

