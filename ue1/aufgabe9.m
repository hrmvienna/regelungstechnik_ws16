
clear all;close all;clc;
% Parameter
Tend = 20;
c = 1; 
m = 1;
k = 2;
A = [0 1; -k/m -c/m];
B = [0;1/m];
C = eye(2); % Damit man die Sprungantwort fuer s und v kriegen
D = 0;
% Anfangswerte
s0 = 1;
v0 = 0.5;
x0 = [s0; v0];

eig_vs = eig(A)
eig_max = max(abs(real(eig_vs)))
eig_min = min(abs(real(eig_vs)))

Steifigkeit = eig_max / eig_min
% nicht steif, da S = 1

m =1; c=1001; k = 1000;
A2 = [0 1; -k/m -c/m];
B2 = [0;1/m];

eig_vs2 = eig(A2)
eig_max2 = max(abs(real(eig_vs)))
eig_min2 = min(abs(real(eig_vs)))

Steifigkeit2 = eig_max2 / eig_min2
% sehr steif, da S = 1000

% Abtastzeit und t-Spanne #1
Ta = 1/1000;
trange = 0:Ta:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig = 500*ones(length(trange),1);

x_euler_steif = int_euler_1_2(A2, B2, u_sig, x0, Ta, Tend);
int_ls = linspace(0, Tend, Tend/Ta + 1);

figure 
subplot(2,1,1);
plot (int_ls, x_euler_steif(1,:), int_ls, u_sig)
legend('euler', 'u sig')
%xlim([-0.5, 2])
ylim([-2, 2])
title(sprintf('s: Ta = %0.4fs', Ta))
grid on

subplot(2,1,2);
plot (int_ls, x_euler_steif(2,:), int_ls, u_sig)
legend('euler', 'u sig')
ylim([-2, 2])
title(sprintf('v: Ta = %0.4fs', Ta))
grid on

Ta = 1/10;
trange = 0:Ta:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig = 500*ones(length(trange),1);
x_euler_steif2 = int_euler_imp(A2, B2, u_sig, x0, Ta, Tend);
int_ls = linspace(0, Tend, Tend/Ta + 1);

figure 
subplot(2,1,1);
plot (int_ls, x_euler_steif2(1,:), int_ls, u_sig)
legend('imp euler', 'u sig')
%xlim([-0.5, 2])
ylim([-2, 2])
title(sprintf('s imp: Ta = %0.4fs', Ta))
grid on

subplot(2,1,2);
plot (int_ls, x_euler_steif2(2,:), int_ls, u_sig)
legend('imp euler', 'u sig')
ylim([-2, 2])
title(sprintf('v imp: Ta = %0.4fs', Ta))
grid on
