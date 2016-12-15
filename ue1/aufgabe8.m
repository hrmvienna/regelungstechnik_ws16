
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

% Abtastzeit und t-Spanne #1
Ta = 0.05; %Anagabe
trange = 0:Ta:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sin = zeros(length(trange),1);

T = 20;
for k = 1:1:length(trange)
    u_sin(k) = sin(10 * pi*(k-1) *Ta / T);
end

x_heun = int_heun(A, B, u_sin, x0, Ta, Tend);
x_rung = int_runge_kutta(A, B, u_sin, x0, Ta, Tend);
x_euler = int_euler_1_2(A, B, u_sin, x0, Ta, Tend);

%exakt
% Exaktes System (sehr fein abgetastet)
sys = ss(A,B,C,D);

lsim(sys, u_sin, trange, x0)

tranges = [0.2, 0.1, 0.05]

for Taa = tranges

% Abtastzeit und t-Spanne #1
trange = 0:Taa:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sin = zeros(length(trange),1);
for k = 1:1:length(trange)
    u_sin(k) = sin(10 * pi*(k-1) *Taa / T);
end

x_heun = int_heun(A, B, u_sin, x0, Taa, Tend);
x_rung = int_runge_kutta(A, B, u_sin, x0, Taa, Tend);
x_euler = int_euler_1_2(A, B, u_sin, x0, Taa, Tend);
    
int_ls = linspace(0, Tend, Tend/Taa + 1);

figure 
subplot(2,1,1);
plot (int_ls, x_heun(1,:), int_ls, x_rung(1,:), int_ls, x_euler(1,:), int_ls, u_sin)
legend('heun', 'runge kutta', 'euler', 'sinus')
title(sprintf('s: Ta = %0.2fs', Taa))
grid on

subplot(2,1,2);
plot (int_ls, x_heun(2,:), int_ls, x_rung(2,:), int_ls, x_euler(2,:), int_ls, u_sin)
legend('heun', 'runge kutta', 'euler', 'sinus')
title(sprintf('v: Ta = %0.2fs', Taa))
grid on

Taa
pause()
end
