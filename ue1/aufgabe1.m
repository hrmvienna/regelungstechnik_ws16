% ===========
% Aufgabe 1.1:
% ===========

clear all   % Es werden alle Variablen gelöscht
close all   % Es werden alle Fenster geschlossen
clc         % Das Command Window wird zurueckgesetzt

%% 
clear all;close all;clc;

% Aufgabe 1.1.1:

% Matrix A und Vektor b
A = [1, 2, 4; 2, 2, 1; 3, 2, 0];
b = [5; 4; 1];

% Zeit
tic
x1 = inv(A)*b
toc
tic
x2 = mldivide(A, b)
toc

% numerischer Fehler
norm(A*x1 - b)
norm(A*x2 - b)

% mldivide ist schneller (0.000148 vs 0.000070) und numerisch genauer

%% 
clear all;close all;clc;

% Aufgabe 1.1.2:

% rect(x) ~~ A* sum(k=1 bis n) 4/pi*(2k-1) * sin((2k-1)x)

A = 10;
x = 0:1:10;
n = 1;
y = 0

for n=1:1:100
    y = y + A*(4/(pi*(2*n -1))) * sin((2*n -1)*x);
    plot(y)
    title(['rec(x) for n = ', num2str(n)])
    grid on
    pause()
    %pause(0.01)
end

%% 
clear all;close all;clc;

% Aufgabe 1.1.3:

% f(x,y,t) = sin(x*pi/10)*sin(y*pi/10)*|sin(t)|
step_size = 0.2;
range = 0:step_size:10;
%times = [0.1, 0.5, 1, 2, 10, 32];
times = 0.1:0.1:10;

[xx,yy] = meshgrid(range,range); % Koordinatenmatrix erzeugen
figure

for t = times
    zz = sin((xx*pi)/10).*sin((yy*pi)/10).*abs(sin(t));

    % Gitterplot
    mesh(xx,yy,zz)
    box on
    grid on
    zlim([0,1])
    
    suptitle(['t = ', num2str(t)]);
    pause(0.1)
end

%%
% ===========
% Aufgabe 1.2:
% ===========

%% 
clear all;close all;clc;

% Aufgabe 1.2.1

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

% Allgemeine Loesung
sa = dsolve('D2s = -(k/m)*s - (c/m)*Ds +1/m')

% partikulaere Loesung
sp = dsolve('D2s = -(k/m)*s - (c/m)*Ds + 1/m','s(0)=0', 'Ds(0)=0')


% Aufgabe 1.2.2
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



%%
% ===========
% Aufgabe 1.3:
% ===========

%% 
clear all;close all;clc;

syms u x1 x2 x3

% System
% x''' +cos(x)^2 * x'' + x' + e^(-x) *u = 0
% x1 = x            ;   x1' = x2
% x2 = x1' = x'     ;   x2' = x3
% x3 = x2' = x''    ;   x3' = -(cos(x1)^2 *x3 + x2 + e(-x1)*u)
x = [x1, x2, x3]';
xx = [x2; x3; -cos(x1)^2 *x3 - x2 + exp(-x1)*u]

% Simulink Parameter
Ta = 0.05;
% Anfangswerte
x1_0 = 0;
x2_0 = 0;
x3_0 = 0;

% Diff.gl.system
% x''' + cos(x)^2 * x'' + x' + e^(-x)*u = 0
% numerisch berechnen und plotten; @myrigid ist eine .m Datei
[T,Y] = ode45(@myrigid,[0 10],[x1_0 x2_0 x3_0]) % timespan [0 10], Anfangsbed.: [0 0 0]

% numerishce Berechnung ploten
figure
plot(T,Y(:,1),'-',T,Y(:,2),'-',T,Y(:,3),'-')
leg_y1 = sprintf('x1');
leg_y2 = sprintf('x2');
leg_y3 = sprintf('x3');
title('Numerische Berechnung von 1.3 a)')
grid on
legend(leg_y1, leg_y2, leg_y3,'Location', 'SouthEast')


% System um die Ruhelage xr = ur = 0 linearisieren
% Am Zettel


%%
% ===========
% Aufgabe 1.4:
% ===========

%% 
clear all;close all;clc;

syms Ue Us Ua Uc1 Uc2 R1 R2 R3 K

% Idealer nichtinvertierender Operationsverstärker
% Allgemein: Ua = (1 + R2/R1) * Ue => Ua = V*Ue
% hier: V = 1 + ((K-1)*R3)/R3 = 1 + K - 1 = K

% Ua = K * Uc2

x = [Uc1; Uc2];
u = Ue;
d = Us;
y = Ua;

