% ===========
% Aufgabe 1.4:
% ===========

%% 
clear all;close all;clc;

%% Aufgabe 1.4.1

syms Ue Us Ua Uc1 Uc2 R1 R2 R3 K
syms Q1 C1 C2 Ur1 Ur2 Ur4 Uc1_ Uc2_
syms is ie ir2 ic1 ic2

% Idealer nichtinvertierender Operationsverstaerker
% Allgemein: Ua = (1 + R2/R1) * Ue => Ua = V*Ue
% hier: V = 1 + ((K-1)*R3)/R3 = 1 + K - 1 = K
% Ua = K * Uc2

x = [Uc1; Uc2];
u = Ue;
d = Us;
y = Ua;

% U = R*i
% Q = C*U
% i_c = dQ/dt = C*dU/dt = C*U'

% C1(Uc1) = dQ1/dUc1
% => dQ1 = C1(Uc1)*dUc1

% Uc{1,2}_ der Unterstrich bezeichnet die Ableitung
is = Ur4/R2;
ie = Ur1/R1;
ic1 = C1*Uc1_;
ic2 = C2*Uc2_;
Ua = K*Uc2;
y = subs(y); % Ua einsetzen

% Knotengleichungen
k1 = ie - ir2 - ic1;
k2 = ic2 - ir2 - is;

% Maschengleichungen
m1 = Ur1 + Ur2 + Uc2 - Ue;
m2 = Ur4 + Uc2 - Us;
m3 = Ur1 + Uc1 + Ua - Ue;
m4 = -Ur2 + Uc1 + Ua - Uc2;

% is aus m2
Ur4 = solve(m2, Ur4);
is = subs(is);

% ie aus m3
Ur1 = solve(m3, Ur1);
ie = subs(ie);

% Loesen nach Ur2_
ir2 = solve (k2, ir2);
Ur2 = R2 * ir2;

m4 = subs(subs(m4));
Uc2_ = simplify(solve(m4, Uc2_));
ic2 = Uc2_ * C2;


% Loesen nach Ur1_
ie = (ic2 - is) + ic1; %k2 in k1 eingesetzt
Ur1 = R1*ie;
m3 = subs(m3);
Uc1_ = simplify(solve(m3, Uc1_));

% X_punkt = f(x, u, d)
Uc1_
Uc2_

%% 1.4.2 - Linearisieren

% Linearisierung durch Anwendung der Taylor-Formel nach Satz 2.5 und 2.6

% Erzeuge die differenzierten Matrizen
A = [diff(Uc1_, Uc1), diff(Uc1_, Uc2); diff(Uc2_, Uc1), diff(Uc2_, Uc2)]
bu = [diff(Uc1_, Ue); diff(Uc2_, Ue)]
bd = [diff(Uc1_, Us); diff(Uc2_, Us)]

cT = [diff(y, Uc1), diff(y, Uc2)]
du = diff(y, Ue)
dd = diff(y, Us)

%% 1.4.3 - Ruhelagen berechnen

syms UeR UsR

Uc1r = solve(Uc1_ == 0, Uc1);
Uc2r = solve(Uc2_ == 0, Uc2);

% in die zweite Geichung setzten wir Uc1r ein und loesen nach Uc2r auf
Uc2r = solve(Uc2 == subs(Uc2r, Uc1, Uc1r), Uc2);

% Uc2r in die erste Gleichung einsetzen und Uc1r berechnen
Uc1r = subs(Uc1r, Uc2, Uc2r);

Uc1r = simplify(subs(Uc1r, [Ue, Us], [UeR, UsR]));
Uc2r = simplify(subs(Uc2r, [Ue, Us], [UeR, UsR]));

% Ruhelagen
x_r = [Uc1r; Uc2r]
yr = subs(y, Uc2, Uc2r)


%% 1.4.4 - Uebertragungsfunktion

syms UC1_ref C1_ref kC1

% Systemmatrix A berechnen: A = [a11 a12 ; a21 a22]
Q1 = (C1_ref + kC1*((Uc1/2) - UC1_ref))*Uc1;
dQ1 = simplify(diff(Q1, Uc1)) % Erste Ableitung von Q1 nach UC1
ddQ1 = simplify(diff(dQ1, Uc1)) % Zweite Ableitung von Q1 nach UC1
C1 = dQ1;
dC1 = ddQ1;

R1      = 625;    % Widerstand in Ohm
R2      = 3500;   % Widerstand in Ohm
K       = 3;      % Spannungsverstaerkung
C1_ref  = 1e-6;   % Referenz-Kapazitaetswert auf Kennlinie von C1 in F
UC1_ref = -10;    % Referenz-Spannungswert auf Kennlinie von C1 in V
kC1     = 800e-9; % Steigung der Kennlinie von C1 in F/V
C2      = 1e-6; 
UeR = 5;
UsR = 5;

A = double(subs(subs(subs(A), [Uc1, Uc2], [Uc1r, Uc2r])));
bu = double(subs(subs(subs(bu), [Uc1, Uc2], [Uc1r, Uc2r])));
bd = double(subs(subs(subs(bd), [Uc1, Uc2], [Uc1r, Uc2r])));
cT = double(subs(subs(subs(cT), [Uc1, Uc2], [Uc1r, Uc2r])));
du = double(subs(subs(subs(du), [Uc1, Uc2], [Uc1r, Uc2r])));
dd = double(subs(subs(subs(dd), [Uc1, Uc2], [Uc1r, Uc2r])));

% Zustandsraumdarstellung 
sys_gl = ss(A,[bu bd],cT,[du dd]);
G_sys = tf(sys_gl);
G = G_sys(1)  % Uebertragungsfunktion von u zu y
Gd = G_sys(2) % Uebertragungsfunktion von d zu y

%% Aufgabe 1.4.5

V = 1.371e06
m = 9.959e05;
% 1 + 2*xi*(s*T) + (s*T)^2 = s^2 + 1600 s + m
% 1 + 1600/m + 1/m s^2 => T^2 = 1/m => T = +/- sqrt(1/m)
% 2*xi*T = 1600/m => xi = 800/sqrt(m)
xi = 800/sqrt(m)
T = sqrt(1/m)


%% Aufgabe 1.4.6
figure
bode(G, Gd)
grid on
legend('G(s)','Gd(s)','Location', 'SouthEast')

% interpretation am Zettel und in Worten
% tl;dr: Die Stoerung wird besser verstaerkt als das Eingangssignal und hat
% eine bessere Phasenreserve