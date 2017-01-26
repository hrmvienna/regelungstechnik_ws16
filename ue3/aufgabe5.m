clear all; close all; clc

%% Aufgabe 3.5

Ta = 0.01;

Rz_komp_den = [1 14 49];
Rz_komp_num = [0.0128407323859812 0.0194498756069027 1];
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
parKompReg.P = 12.770487192571185;
parKompReg.I = 0.435221741157092;
parKompReg.R_kompz = Rz_komp;