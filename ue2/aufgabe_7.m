clear all; close all; clc

%% Aufgabe 2.7.1

% Das Sensorrauschen wirkt sich wegen dem guten Fuehurungsverhalten direkt
% auf die Ausgangsgroesse aus (4.3.1 Ad(5)), Tny = -Try


% Bei kuerzerer Anstiegszeit (tr < 1) und schnellerer Dynamik (mehr Stellgroesse
% wird benoetigt) wird der Verstaerkungsfaktor dementsprechend groesser und
% wirkt sich das Sensorrauschen staerker auf die Stellgroesse aus. Tnu = Tdu

% Bei groesseren Anstiegszeit (tr > 1) wirkt sich das Sensorrauschen nicht
% so stark auf die Stellgroesse aus.


% tr = 0.5, q = 15 -> b_n = -0.8561
% tr = 4, q = 7 -> b_n = -0.0114
% tr = 1, q = 7 -> b_n = -0.1291


%% Parameter

Ta = 50e-03;

%% Linearisiertes System

Gz_lin_num = [0,0.00624342880213793,0.0284109332528076,0.00845608794085389,9.51565947781945e-06];
Gz_lin_den = [1,-2.71439290270362,2.61753640663902,-0.896725281314584,7.17948644118771e-08];
Gz_lin = tf(Gz_lin_num, Gz_lin_den, Ta);

%% PI-Komp Regler aus Aufgabe 2.5
Rz_num = [0.129052836691455 -0.350640525783673 0.338341661652813 -0.115965884301370];
Rz_den = [1 -2.40425531914894 1.89723856948846 -0.492983250339520];
Rz = tf(Rz_num, Rz_den, Ta);

% Regler mit t_r = 4
num = [0.0114235121008542 -0.0302399436766243 0.0285518173838854 -0.00952261914285287];
den = [1 -2.40425531914894 1.89723856948846 -0.492983250339520];
Rz_tr4 = tf(num, den, Ta);

% Regler mit t_r = 0.5, realisierungspol = 15
num = [0.856133395773891 -2.32771511276413 2.24730725704585 -0.770779483086807];
den = [1 -1.90909090909091 1.11570247933884 -0.206611570247934];
Rz_tr05 = tf(num, den, Ta);

%% 

Regler = Rz; % | Rz_tr4 | Rz | Rz_tr05

% Tny: Strecke vom Sensorrauschen zum Ausgang 
Tnu = -Regler/(1+Regler*Gz_lin)

% Tny: Strecke vom Sensorrauschen zum Eingang
Tny = -(Regler*Gz_lin)/(1+Regler*Gz_lin)

[a, b, k_Tnu] = residue(cell2mat(Tnu.Numerator), cell2mat(Tnu.Denominator));
[a, b, k_Tny] = residue(cell2mat(Tny.Numerator), cell2mat(Tny.Denominator));
k_Tnu
k_Tny