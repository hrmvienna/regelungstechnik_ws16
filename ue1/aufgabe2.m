% ===========
% Aufgabe 1.2:
% ===========

%% 
clear all;close all;clc;

%% Aufgabe 1.2.1

%m s(2) + c s(1)  + k s = f

%x1 = s
%x2 = s(1)

%x1. = x2 
%x2. = 1/m * (f - c*x2 - k+x1)

%[x1., x2.] = [[0, 1], [-k/m, -c/m]][x1, x2] + [0, 1/m]*u
%y = [1,0][x1, x2]

syms c  k  m f x1 x2 s0 v0

% Zustandsvektor, Eingangsvektor
x = [x1; x2];
u = f;

% Systemmatrizen
A = [0 1; -k/m -c/m];
B = [0;1/m];
C = [1,0];
D = 0;

xx = A*x + B*u;
y = C*x + D*u;

% Eigenwerte
lambda = eig(A)

% Eigenewrte l1,2 = -(c +/- sqrt(c^2 - 4*k*m)/(2*m)
% c, k, m > 0: Auswirkung auf die Stabilitaet? TODO (z.B. wenn c^2 < 4*k*m,
% dann ist der realteil fix negativ (-c/2m)
lambda
l1 = lambda(1); % = -0.5 + 1.3229 i = -0.5 + sqrt(-7)/2
l2 = lambda(2); % = -0.5 - 1.3229 i = -0.5 - sqrt(-7)/2

% Realteile der Eigenwerte sind negativ, deswegen ist das System stabil.
real(lambda)

% Transitionsmatrix; expm = e^Matrix
t = sym('t');
phi = expm(A * t)

% Einheitssprung Sigma(t) - wir betrachten nur t >= 0
sig = 1;
f = sig;
u = f;

% Allgemeine Loesung, fuer u = Einheitssprung = 1, da t>0 betrachtet wird
sa = dsolve('D2s = -(k/m)*s - (c/m)*Ds +1/m')

% partikulaere Loesung
sp = dsolve('D2s = -(k/m)*s - (c/m)*Ds + 1/m','s(0)=0', 'Ds(0)=0')


%% Aufgabe 1.2.2 und 1.2.3
% Euler-Verfahren Integration

% Parameter
Tend = 20;
c = 1; 
m = 1;
k = 2;
A = [0 1; -k/m -c/m];
B = [0;1/m];
C = eye(2); % Damit man die Sprungantwort fuer s und v kriegen

% Anfangswerte
s0 = 0;
v0 = 0;
x0 = [s0; v0];

% Abtastzeit und t-Spanne #1
Ta = 0.05; %Anagabe
trange = 0:Ta:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig = ones(length(trange),1);
x_euler_1 = int_euler_1_2(A, B, u_sig, x0, Ta, Tend);

% Abtastzeit und t-Spanne #2
Ta2 = 0.02; %Anagabe
trange2 = 0:Ta2:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig2 = ones(length(trange2),1);
x_euler_2 = int_euler_1_2(A, B, u_sig2, x0, Ta2, Tend);

% Abtastzeit und t-Spanne #3
Ta3 = 0.1; %Anagabe
trange3 = 0:Ta3:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig3 = ones(length(trange3),1);
x_euler_3 = int_euler_1_2(A, B, u_sig3, x0, Ta3, Tend);

% Abtastzeit und t-Spanne #4
Ta4 = 0.2; %Anagabe
trange4 = 0:Ta4:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig4 = ones(length(trange4),1);
x_euler_4 = int_euler_1_2(A, B, u_sig4, x0, Ta4, Tend);

% Abtastzeit und t-Spanne #5
Ta5 = 0.4; %Anagabe
trange5 = 0:Ta5:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig5 = ones(length(trange5),1);
x_euler_5 = int_euler_1_2(A, B, u_sig5, x0, Ta5, Tend);

% Struktur zur einfachen speicherung der Daten erstellen
field1 = 'ta';     value1 = {Ta, Ta2, Ta3, Ta4, Ta5};
field2 = 'trange'; value2 = {trange, trange2, trange3, trange4, trange5};
field3 = 'u_sig';  value3 = {u_sig, u_sig2, u_sig3, u_sig4, u_sig5};
field4 = 'x_euler';value4 = {x_euler_1, x_euler_2, x_euler_3, x_euler_4, x_euler_5};

% Datenstruktur
s_eulers = struct(field1,value1,field2,value2,field3,value3,field4,value4);

% Exaktes kontiunierliches System generieren
sys = ss(A,B,C,D);
% UEbertragungsfunktion des kont. Systems
%ksys = tf(sys);
% Sprungantwort fuer das kont. System
%step(ksys, trange(end)) % Graph
% Diskretes System generieren mit Ta = 0.05
dsys = c2d(sys, Ta, 'zoh');
% Sprungantwort fuer das diskrete System
dstep = step(dsys, trange(end));

% Exaktes System (sehr fein abgetastet)
TaEx = 0.001;
trangeEx = 0:TaEx:Tend;
dstep2 = step(c2d(sys, TaEx, 'zoh'), trangeEx);

% Differnzen Plotten
% Wertevectoren fuer s und v in array speichern
field1 = 's';
field2 = 'v';
field3 = 'linespace';
field4 = 'l';
value1 = {dstep2(:,1)};
value2 = {dstep2(:,2)};
value3 = {linspace(0, Tend, Tend/TaEx + 1)};
value4 = {sprintf('Kont. Sys. Ta = %0.3fs', TaEx)};

value1 = [value1, dstep(:,1)];
value2 = [value2, dstep(:,2)];
value3 = [value3, linspace(0, Tend, Tend/Ta + 1)];
value4 = [value4, sprintf('disk. Sys. Ta = %0.2fs', Ta)];

for i=1:length(s_eulers)
    value1 = [value1, s_eulers(i).x_euler(1,:)'];
    value2 = [value2, s_eulers(i).x_euler(2,:)'];
    value3 = [value3, linspace(0, Tend, Tend/(s_eulers(i).ta) + 1)];
    value4 = [value4, sprintf('Euler Ta = %0.2fs', (s_eulers(i).ta))];
end

% Struktur fuer die Plott daten
s_plots_sv = struct(field1, value1, field2, value2, field3, value3, field4, value4);

% Plot fuer s
figure
tt = 0:Tend;
plot(s_plots_sv(1).linespace, s_plots_sv(1).s, 'b',  ...
    s_plots_sv(2).linespace, s_plots_sv(2).s, 'm', ...
    s_plots_sv(3).linespace, s_plots_sv(3).s, 'b', ...
    s_plots_sv(4).linespace, s_plots_sv(4).s, 'g', ...
    s_plots_sv(5).linespace, s_plots_sv(5).s, 'r', ...
    s_plots_sv(6).linespace, s_plots_sv(6).s, 'k', ...
    s_plots_sv(7).linespace, s_plots_sv(7).s, 'c')

%plot(tt, sin(tt))
grid on
xlabel('Zeit t [s]')
ylabel('Weg s')
title('Graph fuer s')
legend(s_plots_sv(1).l, s_plots_sv(2).l, s_plots_sv(3).l, ...
        s_plots_sv(4).l, s_plots_sv(5).l, s_plots_sv(6).l, ...
        s_plots_sv(7).l, 'Location', 'SouthEast')

% plot fuer v
figure
plot(s_plots_sv(1).linespace, s_plots_sv(1).v, 'b',  ...
    s_plots_sv(2).linespace, s_plots_sv(2).v, 'm', ...
    s_plots_sv(3).linespace, s_plots_sv(3).v, 'b', ...
    s_plots_sv(4).linespace, s_plots_sv(4).v, 'g', ...
    s_plots_sv(5).linespace, s_plots_sv(5).v, 'r', ...
    s_plots_sv(6).linespace, s_plots_sv(6).v, 'k', ...
    s_plots_sv(7).linespace, s_plots_sv(7).v, 'c')
grid on
xlabel('Zeit t [s]')
ylabel('Geschweidigkeit v [m/s]')
title('Graph fuer v')
legend(s_plots_sv(1).l, s_plots_sv(2).l, s_plots_sv(3).l, ...
        s_plots_sv(4).l, s_plots_sv(5).l, s_plots_sv(6).l, ...
        s_plots_sv(7).l, 'Location', 'SouthEast')
