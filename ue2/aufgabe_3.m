clear all; close all; clc

%% Aufgabe 2.3.1: Reduziertes System

% Zustandsgroessen
syms i_GSM phi_GSMP w_GSM w_P M_ext
% GSM
syms L_GSM R_GSM k_GSM J_GSM d_cGSM d_vGSM u_GSM
% P
syms J_P d_cP d_vP d_qP
% GSMP
syms c_GSMP d_GSMP phi_GSMP

i_GSM = (1/R_GSM) * u_GSM - (k_GSM/R_GSM) * w_GSM

x_red = [phi_GSMP w_GSM w_P].';
f_red = simplify([w_GSM - w_P; ...
            -(d_cGSM + c_GSMP*phi_GSMP - i_GSM*k_GSM + d_vGSM*w_GSM + d_GSMP*(w_GSM - w_P))/J_GSM; ...
            -(M_ext + d_cP - c_GSMP*phi_GSMP + d_vP*w_P - d_GSMP*(w_GSM - w_P) + d_qP*w_P^2)/J_P]);
d_x_red = f_red
h_red = [0 0 1]*x_red;
y_red = h_red

%% Aufgabe 2.3.2: Ruhelage berechnen und Linearisieren

syms u_GSM_r M_ext_r w_P_r

x_r = solve(subs(d_x_red, [u_GSM M_ext], [u_GSM_r M_ext_r]), x_red);

% Parameterliste
paralist_1 = [u_GSM_r M_ext_r L_GSM    R_GSM k_GSM J_GSM   d_cGSM d_vGSM J_P     d_cP  d_vP   d_qP c_GSMP d_GSMP];
paralist_2 = [5.6     0       1.4e-3   0.46  0.1   12.4e-3 0.152  1.8e-3 32.5e-3 0.169 2.7e-3 1e-4 0.6822 1e-5];

% Ruhelagen mit eingesetzten Werten
r_red_num = [x_r.phi_GSMP(1) x_r.phi_GSMP(2) ; ...
        x_r.w_GSM(1)    x_r.w_GSM(2); ...
        x_r.w_P(1)      x_r.w_P(2)];
r_red_num = double(simplify(subs(r_red_num, paralist_1, paralist_2)));

% Systemmatrizen des linearisierten reduzierten Systems
A = [diff(f_red, phi_GSMP) diff(f_red, w_GSM) diff(f_red, w_P)];
A = simplify(subs(A, w_P, w_P_r))
bu = [diff(f_red, u_GSM)]
bd = [diff(f_red, M_ext)]
ct = [0 0 1];

%% Aufgabe 2.3.3: Ruhelagen vergleichen

% Die Ruhelagen des reduzierten Systems sind gleich dem des vollstaendigen
% Systems.

% Die Naeherung der singulaeren Stroerungstheorie fuer i_GSM stellt genau den
% Zusammenhang zwischen der Stromgleichung und den anderen drei
% Zustandsgleichungen dar. Deshalb ändert sich beim Gleichungssystem nichts, 
% weil die Stromgleichung in den anderen Zustandgleichungen einfließt und 
% somit die selben Ruhelagen besitzt.


% Die Ruhelagen des reduzierten Systems sind identisch mit den Ruhelagen
% des nicht reduzierten Systems, weil das System vorher um die Ruhelagen
% linearisiert wurde und darauf hin reduziert und durch die Näherung der
% singulären Störungstheorie kann die Stromgleichung sehr gut angenehert
% werden (Siehe Simulationsergebnis 2_3_3_i_scope)
