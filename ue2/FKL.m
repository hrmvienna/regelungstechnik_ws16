%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufgabe 2.1
% FKL Reglerentwurf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anforderungen an den offenen und geschlossenen Kreis
tr=3e-3;
ue=5;


%% Aufgabe 2.1.2 PI Regler

T = 1.002e-3;
T_2 = 8.0863e-04;
xi = 0.802;
V = 1.377;
V_G = V;
V_R = 369.8044;

b = [V_G];
a = [(T^2) 2*xi*T 1];
Gs = tf(b, a)

%Gesamtregler und offener Kreis
%R_PI = V_R*(T_2*s +1)/s;
R_PI = V_R*tf([T_2 1], [1 0]);
L_PI = R_PI*Gs;
bode(L_PI,'r');

%% Aufgabe 2.1.3 PID Regler

% Bleibende Regelabweichung
e_inf = 1e-3;

% Zeitkonstante des Intagralterms
TI = 1e-3;

%Gesamtregler und offener Kreis
% R_PID = 
% L_PID = 
% bode(L_PID,'g');

%% Aufgabe 2.1.4
% Kontrolle des Verhaltens des geschlossenen Kreise

T = 0.05;   % Simulationsdauer
figure;
subplot(3,1,1);hold on;grid on;
title('Sprungantwort Führungsübertragungsfunktion')

subplot(3,1,2);hold on;grid on;
title('Sprungantwort Störübertragungsfunktion')

subplot(3,1,3);hold on;grid on;
title('Rampenantwort Störübertragungsfunktion')