clear all; close all; clc

%% Aufgabe 3.5
% save('regler14.mat', 'parKompReg', 'parZR')

Ta = 0.01;

% aus Aufgabe 3.2:
A = double([0, 1, -1; -3411/62, -54163/28520, 1/1240; ...
    6822/325, 1/3250, 54117/74750 - (7*1130385^(1/2))/7475]);
bu = double([0; 12500/713; 0]);
bd = double([0; 0; -400/13]);
ct = [0 0 1]; d = 0;
sys = ss(A,bu, ct, d);
dsys = c2d(sys, Ta);

Rz_komp_den = [1 -1.86473429951691 0.869308501948704];
Rz_komp_num = [0.0124660457258177 -0.0246598980273402 0.0122872033715591];
Rz_komp = tf(Rz_komp_num, Rz_komp_den, Ta);

% PI-Zustandsregler aus Aufgabe 3.3

kI = 0.0065;
kP = 0.1488;
kx = [-0.6050   -0.7335    0.3827];

% Parameter fuer PI Regler als Matlab Function
parZR.c_z = dsys.c;
parZR.kI = kI;
parZR.kp = kP;
parZR.kx = kx;
parZR.Ta = Ta;

% Parameter Kompensationsregler
parKompReg.P = 12.432087412158705;
parKompReg.I = 0.090619696734693;
parKompReg.R_kompz = Rz_komp;