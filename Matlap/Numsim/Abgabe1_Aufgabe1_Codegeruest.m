% Bestimme Naeherungsloesung fuer ein Feder-Daempfer-Masse-System:
% m*y''(t) + D*y'(t) + k*y(t) = F_0*cos(omega*t-Delta)
% mit y(t=0)=y_0, y'(t=0)=ydot_0 und m,k,D = const 

%% HAUPTFUNKTION
% =========================================================================
function AWP_Loeser()

format shortE
close all

% Modellparameter
m     = % Masse
k     = % Federhaerte
D     = % Daempfung
omega = % Anregungskreisfrequenz
Delta = % Phasenwinkel
F_0   = % Kraftamplitude

% Simulationsparameter
y_0   = % Anfangswerte
t_0   = % Startzeitpunkt
t_end = % Endzeitpunkt

% Diskretisierungsparameter
N       = % Anzahl Intervallschritte
h       = % Schrittweite
t       = % Laufvariable
yh      = % Loesungsvariable
n_start = % Startindex Zeititeration

% Loesungsvektor anlegen
y_numerisch      = zeros(3,N+1);    % erste Zeile = Stuetzstellen t, zweite Zeile = Loesungswerte y1, dritte Zeile = Loesungswerte y2
y_numerisch(:,1) = [t_0, y_0(1), y_0(2)];

% Zeit-Iterationsschleife .............................................
for n = n_start:N+1
  % Berechnung Loesungsinkrement
  dy = RK4(...);

  % Lauf- und Loesungsvariable inkrementieren
  % ...

  % Werte in Loesungsvektor speichern
  y_numerisch(:,n) = [t;yh'];
end % .................................................................

% quadratischen Fehler berechnen und ausgeben
% ...
    
% analytische und numerische Loesung plotten
% ...

end % Funktion AWP_Loeser()


%% UNTERFUNKTIONEN
% =========================================================================

% rechte Seite (RHS) der DGL ----------------------------------------------
function [RHS] = f(...)
  % ...
end

% Erregungskraft ----------------------------------------------------------
function [F] = anregung(...)
    % ...
end

% Loesungsinkrement dy des klassischen Runge-Kutta-Verfahrens -------------
function [dy] = RK4(...)
  % ...
end

% analytische Loesung y(t) fuer erzwungene Schwingung ---------------------
function [output] = y1_erzwungen_exakt(...)
  % ...
end