clear all; close all; clc

%% Aufgabe 3.3: Zeitdiskreter PI-Zustandsregler

% Lin. red. Modell aus Aufgabe 2.3
A = double([0, 1, -1; -3411/62, -54163/28520, 1/1240; ...
    6822/325, 1/3250, 54117/74750 - (7*1130385^(1/2))/7475]);
bu = double([0; 12500/713; 0]);
bd = double([0; 0; -400/13]);
ct = [0 0 1]; d = 0;
sys = ss(A,bu, ct, d);

% Abtastsystem
Ta = 10e-03; % 10ms
dsys = c2d(sys, Ta);
phi = dsys.A;
gamma = dsys.B;

%% PI-Zustandsregler - Skript Kapitel 8.2, (8.52)
AI = [phi,[0,0,0]'; -ct, 1];
bI = [gamma; 0];

% Gewuenschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
lambda0 = -4; PI = [lambda0, lambda0, lambda0, lambda0];
PdI = exp(PI*Ta); % gew. Pole des geschl. Kreises für das Abtastsystem

% Polvorgabe
ke = -acker(AI,bI,PdI);

%Test
eig(AI+bI*ke)

%Berechnung der Regleranteile gemaess Skriptum
kI = ke(4);
kP = 1/(ct*inv(eye(3)-phi)*gamma);
kx = kP*ct+ke(1:3);

%% Parameter für PI Regler als Matlab Function
parZR.Phi      = dsys.a;
parZR.Gamma    = dsys.b;
parZR.C        = dsys.c;
parZR.D        = dsys.d;

parZR.kI = kI;
parZR.kP = kP;
parZR.kx = kx;


%% Parameter für NLIN - System aus Aufgabe4

u_GSM_r = 5.6;
M_ext_r = 0;

% Parameter
L_GSM   = 1.4e-3;
R_GSM   = 0.46;
k_GSM   = 0.1;
J_GSM   = 12.4e-03;
d_cGSM  = 0.152;
d_vGSM  = 1.8e-03;
J_P     = 32.5e-03;
d_cP    = 0.169;
d_vP    = 2.7e-03;
d_qP    = 1e-04;
c_GSMP  = 0.6822;
d_GSMP  = 1e-05;

% TODO Angaben anpassen, falls sich was in Aufgabe 2.2/2.3 aendert
% Ruhelagen
i_GSM_r = (u_GSM_r + (k_GSM*(R_GSM*d_vP + R_GSM*d_vGSM - ... 
    (k_GSM^4 + R_GSM^2*d_vP^2 + R_GSM^2*d_vGSM^2 - ...
    4*M_ext_r*R_GSM^2*d_qP - 4*R_GSM^2*d_cP*d_qP - ...
    4*R_GSM^2*d_cGSM*d_qP + 2*R_GSM^2*d_vP*d_vGSM + ...
    2*R_GSM*d_vP*k_GSM^2 + 2*R_GSM*d_vGSM*k_GSM^2 + ...
    4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + k_GSM^2))/(2*R_GSM*d_qP))/R_GSM;

phi_GSMP_r = ((k_GSM^2 + R_GSM*d_vGSM)*(R_GSM*d_vP + R_GSM*d_vGSM - ...
    (k_GSM^4 + R_GSM^2*d_vP^2 + R_GSM^2*d_vGSM^2 - ...
    4*M_ext_r*R_GSM^2*d_qP - 4*R_GSM^2*d_cP*d_qP - ...
    4*R_GSM^2*d_cGSM*d_qP + 2*R_GSM^2*d_vP*d_vGSM + ...
    2*R_GSM*d_vP*k_GSM^2 + 2*R_GSM*d_vGSM*k_GSM^2 + ...
    4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + k_GSM^2))/(2*R_GSM^2*c_GSMP*d_qP)...
    - (R_GSM*d_cGSM - k_GSM*u_GSM_r)/(R_GSM*c_GSMP);

w_GSM_r = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

w_P_r = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

w_P_r_red = -(R_GSM*d_vP + R_GSM*d_vGSM - (k_GSM^4 + R_GSM^2*d_vP^2 + ...
    R_GSM^2*d_vGSM^2 - 4*M_ext_r*R_GSM^2*d_qP - ...
    4*R_GSM^2*d_cP*d_qP - 4*R_GSM^2*d_cGSM*d_qP + ...
    2*R_GSM^2*d_vP*d_vGSM + 2*R_GSM*d_vP*k_GSM^2 + ...
    2*R_GSM*d_vGSM*k_GSM^2 + 4*R_GSM*d_qP*k_GSM*u_GSM_r)^(1/2) + ...
    k_GSM^2)/(2*R_GSM*d_qP);

% Anfangsbedingungen fuers Nicht-Lineare System
i_GSM_0 = i_GSM_r;
phi_GSMP_0 = phi_GSMP_r;
w_GSM_0 = w_GSM_r;
w_P_0 = w_P_r;

x_r = [i_GSM_r; phi_GSMP_r; w_GSM_r; w_P_r];