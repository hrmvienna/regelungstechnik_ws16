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

x1r = 0;
x2r = 0;
x3r = 0;
ur = 0;

% Taylorform 2. Ordnung, laut Satz 2.5 u. 2.6
A = [diff(xx, x1) diff(xx, x2) diff(xx, x3)]
b = [diff(xx, u)]
% Ruhelagen einsetzen
A = subs (A, [x1, x2, x3, u], [x1r, x2r, x3r, ur])
b = subs (b, [x1], [x1r])
% Realteile der Eigenwerte von A, zur bestimmung der globalen
% asymptotischen Stabilitaet
real(eig(A))