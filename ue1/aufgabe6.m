% ===========
% Aufgabe 1.6:
% ===========

%% 
clear all;close all;clc;

%% Aufgabe 1.6.1 - Zustandsraumdarstellung des Prozessors

syms dc dck dp dpk c T1 s Ta Tau
syms A b ct d G_s Phi_s Gamma_s Phi_p Gamma_p

% G(s) = dc_s / dp_s = (c)/(s*T1 + 1) = (c/t1)*(1/(s + 1/T1))
% laut Satz 3.5 (Seite 57):
% G(s) = c'*(sE -A)^(-1)*b + d

% Phi_p = (sE -A)^(-1) = 1/(s + 1/T1)
% Koeffizientenvergleich, Bzw. ueber die 2. std. Normalform
% => s - A = s + 1/t1 => A = -1/T1
% => c' = 1 , b = c/T1 , d = 0

A = -1/T1;
b = c/T1;
ct = 1;
d = 0;
Phi_s = (s-A)^(-1);

G_s = simplify(ct*Phi_s*b + d)


% Kontinuierliches System, 2. std. Normalform des Prozessors
% mit u = dp und y = dc
dc_ = A*dc + b*dp
y   = ct*dc + d*dp

% Abtastsystem dazu 
% wie auf Seite 145 auf ueber Satz 2.4 hergeleitet, Gleichung 6.19a,b
% Gamma_p = int (0, Ta) [exp(A*t) dt] = Ta(1 - e^(-Ta/T1))

Phi_p = exp(A*Ta);
Gamma_p = int(exp(A*Tau), Tau, 0, Ta)*b


dck_1 = Phi_p*dck + Gamma_p*dpk
y_p = ct*dck

%% Aufgabe 1.6.2 - Gesamtsystem

% uk = din - dout = dpk
% dpk = a*Sk
% yk = S
syms a Sk uk din dout
% Zustandsgroessen fuer das Gesamtsystem: x_g = [Sk, dck]
% Sk+1 = aktueller Inhalt + ("Zufluesse" - "Abfluesse") * Abtastzeit [Da 
% die Zufluesse und Abluesse als Raten angegeben sind]
% Sk+1 = Sk + Ta*uk + Ta*dck - Ta*dpk
% Sk+1 = Sk + Ta*(uk + dck - a*Sk);
% dck+1 = Phi_p * xk + Gamma_p*dpk
dpk = a*Sk;

% Gesamtsystemgleichungen
Sk_1  = (1 - a*Ta)*Sk + Ta*(uk + dck)
dck_1 = Phi_p*dck + Gamma_p*dpk

% Matrizen und Zustandsvektor der Zustandsraumdarstellung fuer 
% die Zustaende Sk und dck; _g steht fuer Gesamt
x_g= [Sk; dck]

A_g = [(1-Ta*a) Ta; (1-exp(-Ta/T1))*a*c exp(-Ta/T1)]
B_g = [Ta ; 0];
C_g = [1 0];
D_g = 0;

%% Aufgabe 1.6.3

% Parameter
Ta = 10^(-6); % 1 mus
T1 = 10^(-3); % 1 ms
c = 0.5;
a = 0.4; %s^-1

Phi_p_v = double(subs(Phi_p))
Gamma_p_v = double(subs(Gamma_p))
ct_v = double(subs(ct))
A_g_v = double(subs(A_g))
B_g_v = double(subs(B_g))

% Diskrete Zustandsdarstellung fuer den Processor
sys_p = ss(Phi_p_v, Gamma_p_v, ct_v, d, Ta)
%Gz = tf(sys_p);

% Prozessor G(s) mit system toolbox in den Zustandsraum transformiert
%Gs_p = tf(c, [T1, 1]);
%Sys_p = ss(Gs_p);
%Dsys_p = c2d(Sys_p, Ta)
% Dsys_p == sys_p => korrekt

%% Aufgabe 1.6.4

% Ruhelage Sr fuer ur = 1
% x_ = A_g * x_gr + B*ur = 0
% => x_gr = -A^-1 * B*ur
ur = 1;
x_gr = - mldivide(A_g, B_g)*ur

% TODO mit endwerdsatz berechnen ?

%% Aufgabe 1.6.5

% TODO