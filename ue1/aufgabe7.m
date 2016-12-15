
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
x0 = [0; 0];

% Abtastzeit und t-Spanne #1
Ta = 0.05; %Angabe
trange = 0:Ta:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig = ones(length(trange),1);

x_heun = int_heun(A, B, u_sig, x0, Ta, Tend);
x_rung = int_runge_kutta(A, B, u_sig, x0, Ta, Tend);
%x_euler = int_euler_1_2(A, B, u_sig, x0, Ta, Tend);

%exakt
% Exaktes System (sehr fein abgetastet)
% bei ss kann man keinen Startwert angeben!
sys = ss(A,B,C,D);
TaEx = 0.01;
trangeEx = 0:TaEx:Tend;
dstep2 = step(c2d(sys, TaEx, 'zoh'), trangeEx);

tranges = [0.2, 0.1, 0.05]

initial(sys,x0)

for Taa = tranges

% Abtastzeit und t-Spanne #1
trange = 0:Taa:Tend;
% Einheitssprung Sigma(t) als Eingangssignal:
u_sig = ones(length(trange),1);

x_heun = int_heun(A, B, u_sig, x0, Taa, Tend);
x_rung = int_runge_kutta(A, B, u_sig, x0, Taa, Tend);
%x_euler = int_euler_1_2(A, B, u_sig, x0, Taa, Tend);
    
int_ls = linspace(0, Tend, Tend/Taa + 1);
ext_ls = linspace(0, Tend, Tend/TaEx + 1);

figure 
plot (int_ls, x_heun(1,:), ext_ls, dstep2(:,1)')
legend('heun', 'exakt')
title(sprintf('Ta = %0.2fs', Taa))
grid on

figure
plot(int_ls, x_rung(1,:), ext_ls, dstep2(:,1)')
legend('rung','exakt')
title(sprintf('Ta = %0.2fs', Taa))
grid on
Taa
pause()
end

Taa = 0.05
x_heun = int_heun(A, B, u_sig, x0, Taa, Tend);
%x_rung = int_runge_kutta(A, B, u_sig, x0, Taa, Tend);
x_euler = int_euler_1_2(A, B, u_sig, x0, Taa, Tend);
int_ls = linspace(0, Tend, Tend/Taa + 1);


figure
plot(int_ls, x_heun(1,:), int_ls, x_euler(1,:))
legend('heu','euler')
title(sprintf('Ta = %0.2fs', Taa))
grid on
Taa

