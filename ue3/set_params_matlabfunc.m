%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Institut für Automatisierungstechnik 
% Gruppe für komplexe dynamische Systeme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all 
clc

% Kontinuierliches Testsystem 
A = [0 1; -1 -1];
B = [0;1];
C = [1;0];
D = 0;
sys = ss(A,B,C',D);

% Diskretes Testsystem 
Ta = 0.1;
sysd = c2d(sys,Ta);

% Parameter für Matlab Function
parSys.Phi      = sysd.a;
parSys.Gamma    = sysd.b;
parSys.C        = sysd.c;
parSys.D        = sysd.d;

% Anfangszustand
x0 = [0.2; 0.4];
