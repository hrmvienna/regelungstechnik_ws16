% ===========
% Aufgabe 1.6:
% ===========

%% 
clear all;close all;clc;

%% Aufgabe 1.6.1 - Zustandsraumdarstellung des Prozessors

syms dc dck dp dpk c T1 s Ta
syms A b ct d G_s Phi_s Gamma_s Phi_p Gamma_p

% G(s) = dc_s / dp_s = (c)/(s*T1 + 1) = (c/t1)*(1/(s + 1/T1))
% laut Satz 3.5 (Seite 57):
% G(s) = c'*(sE -A)^(-1)*b + d

% Phi_p = (sE -A)^(-1) = 1/(s + 1/T1)
% Koeffizientenvergleich
% => s - A = s + 1/t1 => A = -1/T1
% => c' = c/T1 , b = 1 , d = 0

A = -1/T1;
b = c;
ct = 1;
d = 0;
Phi_s = (s-A)^(-1);

G_s = simplify(ct*Phi_s*b + d)

% Kontinuierliches System
% mit u = dp und y = dc
dc_ = A*dc + b*dp
y   = ct*dc + d*dp

% Abtastsystem dazu 
% wie auf Seite 145 auf Satz 2.4 angewandt
% Gamma_p = int (0, Ta) [exp(A*t) dt] = Ta(1 - e^(-Ta/T1))

Phi_p = exp(A*Ta);
Phi_p_ = diff (Phi_p, Ta);
Gamma_p = (1 - Phi_p)*b;

dck_ = Phi_p*dck + Gamma_p*dpk
y_p = ct*dck

%simplify(subs(dck_, [Ta, T1, c], [1*10^(-9), 1*10^(-3), 0.5]));

%% Aufgabe 1.6.2

% uk = din - dout = dpk
% dpk = a*Sk
% yk = S
syms a Sk uk din dout
% Zustandsgrößen für das Gesamtsystem: xg = [Sk, dck]
% Sk+1 = aktueller Inhalt + ("Zufluesse" - "Abfluesse") * Abtastzeit [Da 
% die Zufluesse und Abluesse als Raten angegeben sind]
% Sk+1 = Sk + Ta*uk + Ta*dck - Ta*a*Sk
% Sk+1 = Sk + Ta*(uk + dck - dpk);
% dck+1 = ct*xk+1
% dck+1 = ct*(Phi_p * xk + Gamma_p*dpk)
dpk = a*Sk;

% Gesamtsystemgleichungen
Sk_  = (1 - a*Ta)*Sk + Ta*(uk + dck)
dck_ = Phi_p*dck + Gamma_p*dpk

%% Aufgabe 1.6.3
set_param('Aufgabe_1_6','AlgebraicLoopSolver','LineSearch')
% Parameter
Ta = 10^(-6); % 1 mus
T1 = 10^(-2); % 1 ms
c = 0.5;
a = 0.4; %s^-1

Phi_p_v = double(subs(Phi_p))
Gamma_p_v = double(subs(Gamma_p))
b_v = double(subs(b))
A_v = double(subs(A))
ct_v = double(subs(ct))

sys_p = ss(Phi_p_v, Gamma_p_v, ct_v, d, Ta)
Gz = tf(sys_p)

% Prozessor G(s) mit system toolbox in den Zustandsraum transformiert
Gs_p = tf(c, [T1, 1])
ss(Gs_p)

%% Aufgabe 1.6.4

% TODO

%% Aufgabe 1.6.5

% TODO