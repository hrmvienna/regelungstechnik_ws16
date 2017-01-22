% 
% Automatisierung: Kapitel 8
% Zum Zustandsregler- und Zustandsbeobachterentwurf in Matlab
% F�r weitere Details konsultieren Sie bitte die help-Funktion von Matlab
%
% Zeitkontinuierliches Modell f�r die SIMULINK-Simulation
% Stand: WS 2016/2017
clear all;
clc;
%% Simulationsmodell
Ac  = [0,1,0,0;-1,-1,1,1;0,0,0,1;0.1,0.1,-0.2,-0.2];
bc  = [0,-1,0,0]';
bvc = [0,0,0,-0.1]';
Bc  = [bc,bvc];
%Messung aller Zustandsgr��en
Cc  = eye(4);
Dc  = [0,0;0,0;0,0;0,0];

%% Entwurfsmodell
Ac;
bc;
cc = [0,0,1,0];
dc = 0;
sysC = ss(Ac,bc,cc,dc);

%Abtastzeit Ta 
Ta = 2;

%Abtastmodell
sysD = c2d(sysC,Ta,'zoh');
Ad = sysD.a;
bd = sysD.b;
%% Reglerentwurf durch Polvorgabe im Zustandsraum
%(File: Automatisierung_Kapitel8_zustandsregler.mdl)
%
%Eigenwerte der Dynamikmatrix im Zeitkontinuierlichen
eig(Ac)

%Gew�nschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
pc = [-0.5+0.5*1i,-0.5-0.5*1i,-1,-2];

%Gew�nschte Pole des geschlossenen Kreises f�r das Abtastsystem
pd = exp(pc*Ta);

%Polvorgabe mit der Ackermannformel (ACHTUNG Ackermannformel in Matlab liefert ein anderes 
%Vorzeichen als im Skriptum!!!)
k = -acker(Ad,bd,pd);

%Test
eig(Ad+bd*k)

%Vorfaktor
g = 1/((cc+dc*k)*inv(eye(4)-Ad-bd*k)*bd+dc);

%% Reglerentwurf durch Polvorgabe im Zustandsraum mit Integralanteil
%(File: Automatisierung_Kapitel8_zustandsregler_Integralanteil.mdl)
%
%Erweiterte Dynamikmatrix und erweiterter Eingangsvektor
AdI = [Ad,[0,0,0,0]';-cc,1];
bdI = [bd;0];

%Gew�nschte Pole des geschlossenen Kreises im Zeitkontinuierlichen
pcI = [-0.5+0.5*1i,-0.5-0.5*1i,-1,-2,-3];

%Gew�nschte Pole des geschlossenen Kreises f�r das Abtastsystem
pdI = exp(pcI*Ta);

%Polvorgabe mit der Ackermannformel (ACHTUNG Ackermannformel in Matlab liefert ein anderes 
%Vorzeichen als im Sktiptum!!!)
ke = -acker(AdI,bdI,pdI);

%Test
eig(AdI+bdI*ke)

%Berechnung der Regleranteile gem�� Skriptum
kI = ke(5);
kP = 1/(cc*inv(eye(4)-Ad)*bd);
kx = kP*cc+ke(1:4);

%% Beobachterentwurf durch Polvorgabe im Zustandsraum
%(File: Automatisierung_Kapitel8_zustandsbeobachter.mdl)
%

%Anfangswert zum Test des trivialen Beobachters
x0 = [10,0,0,0];

%Trivialer Beobachter
AdT = Ad;
bdT = bd;
CdT = eye(4);
ddT = [0,0,0,0]';

%Gew�nschte Pole der Fehlerdynamik im Zeitkontinuierlichen
pcbeo = [-3-3*i,-3+3*i,-1+i, -1-i];

%Gew�nschte Pole des geschlossenen Kreises f�r das Abtastsystem
pdbeo = exp(pcbeo*Ta);

%Entwurf des vollst�ndigen Luenberger Beobachters durch Polvorgabe mit der Ackermannformel 
%(ACHTUNG Ackermannformel in Matlab liefert ein anderes Vorzeichen als im Sktiptum!!!)
kbeo = -acker(Ad',cc',pdbeo);

%Vollst�ndiger Luenberger Beobachter
AdL = Ad+kbeo'*cc;
BdL = [bd,-kbeo'];
CdL = eye(4);
DdL = [0,0;0,0;0,0;0,0];
