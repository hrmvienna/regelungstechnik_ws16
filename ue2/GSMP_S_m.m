function GSMP_S_m(block)

% Simulationsmodell fuer die Gleichstrommaschine (GSM) mit Propeller (P)
%
% -------------------------------------------------------------------------
%
% Beschreibung: TODO
%
% -------------------------------------------------------------------------
%
% Eingaenge:    u1(1) ... u_GSM   Eingangsspannung [V]
%               u2(1) ... M_ext   externes Moment [Nm]
%
% Zustaende:    x(1)  ... i_GSM    Strom an der Gleichstrommaschine [A]
%               x(2)  ... phi_GSMP Drehwinkeldifferenz der GSM und des
%                                  Propellers [rad]
%               x(3)  ... w_GSM    Winkelgeschwindigkeit der GSM [rad/s]
%               x(4)  ... w_P      Winkelgeschwindigkeit des P [rad/s]
%
% Ausgaenge:    y1(1) ... i_GSM    Strom an der Gleichstrommaschine [A]
%               y1(2) ... phi_GSMP Drehwinkeldifferenz der GSM und des
%                                  Propellers [rad]
%               y1(3) ... w_GSM    Winkelgeschwindigkeit der GSM [rad/s]
%               y1(4) ... w_P      Winkelgeschwindigkeit des P [rad/s]
%
%               y2(1) ... M_GSM    elektrische Moment der GSM [Nm]
%               y2(2) ... M_kopp   Kopplungsmoment [Nm]
%
% Parameter: 
%           p(1)... L_GSM       TODO
%           p(2)... R_GSM       TODO
%           p(3)... k_GSM       TODO
%           p(4)... J_GSM       TODO
%           p(5)... d_cGSM      TODO
%           p(6)... d_vGSM      TODO
%           p(7)... J_P         TODO
%           p(8)... d_cP        TODO
%           p(9)... d_vP        TODO
%           p(10)... d_qP       TODO
%           p(11)... c_GSMP     TODO
%           p(12)... d_GSMP     TODO
%           p(13)... i_GSM_0    Anfangszustand i_GSM
%           p(14)... phi_GSMP_0 Anfangszustand phi_GSMP
%           p(15)... w_GSM_0    Anfangszustand w_GSM
%           p(16)... w_P        Anfangszustand w_P
%
% -------------------------------------------------------------------------
% Abtastzeit (sample time): zeitkontinuierlich (continuous)
% -------------------------------------------------------------------------


% Die Funktion setup (s.u.) dient der Initialiserung des Matlab Objektes
% (block). Im Objekt block sind alle fuer die Simulation in Simulink
% notwendigen Eigenschaften (Eingaenge, Zustaende, Ausgaenge, Parameter,
% usw.) des dynamischen Systems (math. Modell) zusammengefasst.
setup(block);

% -------------------------------------------------------------------------
% Initialisierung des Simulationsobjektes block
% -------------------------------------------------------------------------

function setup(block)
  
  % Anzahl der Eingangs- und Ausgangsports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 2;
  
  % Anzahl der zeitkontinuierlichen Zustaende
  block.NumContStates = 4;

  % Anzahl der Parameter
  block.NumDialogPrms = 16;
  
  % Dimensionen der Eingangsports
  % Flag DirectFeedthrough kennzeichnet, ob ein Eingang direkt an einem
  % Ausgang auftritt, d.h. y=f(u)
  block.InputPort(1).Dimensions        = 2;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).DirectFeedthrough = false;

  % Dimensionen der Ausgangsports  
  block.OutputPort(1).Dimensions       = 4;
  block.OutputPort(1).SamplingMode = 'Sample';
  block.OutputPort(2).Dimensions       = 2;
  block.OutputPort(2).SamplingMode = 'Sample';
  
  
  % Einstellen der Abtastzeit: [0 0] wird verwendet fuer die
  % zeitkontinuierliche Simulation.
  block.SampleTimes = [0 0];
  
  % ------------------------------------------------
  % NICHT VERAENDERN
  % ------------------------------------------------
  % 
  % Registrieren der einzelnen Methoden
  % Hier: InitializeConditions ... Initialisierung
  %       Outputs ...       Berechnung der Ausgaenge
  %       Derivatives ...   Berechnung der Zustaende
  %       Terminate ...     Konsistentes Beenden der Simulation

  block.RegBlockMethod('InitializeConditions',    @InitConditions); 
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivatives);  
  block.RegBlockMethod('Terminate',               @Terminate);


% -------------------------------------------------------------------------
% Setzen der Anfangsbedingungen der Zustaende
% -------------------------------------------------------------------------

function InitConditions(block)
  
  % Einlesen der Parameter des Systems
  i_GSM_0    = block.DialogPrm(13).Data;
  phi_GSMP_0 = block.DialogPrm(14).Data;
  w_GSM_0    = block.DialogPrm(15).Data;
  w_P_0      = block.DialogPrm(16).Data;
  
  % Eingabe der Anfangsbedingungen
  x0(1) = i_GSM_0;
  x0(2) = phi_GSMP_0;
  x0(3) = w_GSM_0;
  x0(4) = w_P_0;
  
  % Schreiben auf Objekt block (NICHT VERAENDERN)
  block.ContStates.Data = x0;


% -------------------------------------------------------------------------
% Berechnen der Ausgaenge
% -------------------------------------------------------------------------

function Output(block)

  % Einlesen der Parameter des Systems - TODO: unnoetige entfernen
  L_GSM   = block.DialogPrm(1).Data;
  R_GSM   = block.DialogPrm(2).Data;
  k_GSM   = block.DialogPrm(3).Data;
  J_GSM   = block.DialogPrm(4).Data;
  d_cGSM  = block.DialogPrm(5).Data;
  d_vGSM  = block.DialogPrm(6).Data;
  J_P     = block.DialogPrm(7).Data;
  d_cP    = block.DialogPrm(8).Data;
  d_vP    = block.DialogPrm(9).Data;
  d_qP    = block.DialogPrm(10).Data;
  c_GSMP  = block.DialogPrm(11).Data;
  d_GSMP  = block.DialogPrm(12).Data;
   
  % Shortcut fuer den Zustand
  x = block.ContStates.Data;
  i_GSM    = x(1);
  phi_GSMP = x(2);
  w_GSM    = x(3);
  w_P      = x(4);

  % Berechnung der Ausgaenge
  % Port 1:
  y1(1) = i_GSM;
  y1(2) = phi_GSMP;
  y1(3) = w_GSM;
  y1(4) = w_P;
  
  % Port 2: [M_GSM M_kopp]
  M_GSM = k_GSM * i_GSM;
  M_kopp = (w_GSM - w_P) * d_GSMP + (phi_GSMP)*c_GSMP;
  
  y2(1) = M_GSM;
  y2(2) = M_kopp;
  
  % Schreiben auf Objekt block
  block.OutputPort(1).Data = y1;
  block.OutputPort(2).Data = y2;
  

% -------------------------------------------------------------------------
% Berechnen der Zustaende
% -------------------------------------------------------------------------

function Derivatives(block)

  % Einlesen der Parameter des Systems - TODO: unnoetige entfernen
  L_GSM   = block.DialogPrm(1).Data;
  R_GSM   = block.DialogPrm(2).Data;
  k_GSM   = block.DialogPrm(3).Data;
  J_GSM   = block.DialogPrm(4).Data;
  d_cGSM  = block.DialogPrm(5).Data;
  d_vGSM  = block.DialogPrm(6).Data;
  J_P     = block.DialogPrm(7).Data;
  d_cP    = block.DialogPrm(8).Data;
  d_vP    = block.DialogPrm(9).Data;
  d_qP    = block.DialogPrm(10).Data;
  c_GSMP  = block.DialogPrm(11).Data;
  d_GSMP  = block.DialogPrm(12).Data;
  
  % Shortcut fuer den Eingang
  u = block.InputPort(1).Data;
  u_GSM = u(1);
  M_ext = u(2);
  
  % Shortcut fuer die Zustaende
  x = block.ContStates.Data;
  i_GSM    = x(1);
  phi_GSMP = x(2);
  w_GSM    = x(3);
  w_P      = x(4);
  
  M_GSM = k_GSM * i_GSM;
  M_rGSM = d_cGSM + d_vGSM * w_GSM;
  
  % Berechnen der Zeitableitungen der Zustaende
  dx(1) = (1/L_GSM)*(u_GSM - R_GSM*i_GSM - k_GSM*w_GSM);
  dx(2) = w_GSM - w_P;
  dx(3) = (1/J_GSM)*(M_GSM - M_rGSM - M_kopp);
  dx(4) = (1/J_P)*(M_kopp - M_P);
  
  % Schreiben auf Objekt block
  block.Derivatives.Data = dx;


% -------------------------------------------------------------------------
% Operationen am Ende der Simulation
% -------------------------------------------------------------------------

% Die function Terminate wird hier nicht verwendet,
% muss aber vorhanden sein!
function Terminate(block)

