%% HAUPTFUNKTION
% Bestimme Naeherungsloesung fuer ein Feder-Daempfer-Masse-System:
% m*y''(t) + D*y'(t) + k*y(t) = 0   , y(t=0)=y_0, y'(t=0)=0
% mit m,k,D = const 
% =========================================================================
function AWP_Loeser()

% format shortE
 close all

% Modell-Parameter


% Integrationsverfahren:
% - Heun:                     'HN'  (Aufgabe 1)
% - Adams-Bashforth:          'AB'  (Aufgabe 2)
% - Adams-Bashforth-Moulton:  'ABM' (Aufgabe 3)
integrator = 'HN';

% Anfangswert, Start-/Endzeitpunkt
x0 = 0;
x_end = 1;

t0 = 0;
t_end = 10;
D = 0.1;
k = 1;
N = 40;
h = (t_end - t0)/N;
m = 1 ; 
y0 = 1;

% Anzahl Intervallschritte, Schrittweite

% Loesungsvektor anlegen
y_numerisch = zeros(3,N+1);    % erste Zeile = Stuetzstellen t, zweite Zeile = Loesungswerte y1, dritte Zeile = Loesungswerte y2

% Anfangswert in Loesungsvektor speichern
y_numerisch(:,1) =  [t0;y0;0];
    
% Initialisierung Lauf- und Loesungsvariab  le
switch integrator
  case 'HN'
    % Heun-Verfahren (Aufgabe 1): erste Stuetzstelle = Anfangswert
    t  = t0;
    yh =  [y0 0];
    n_start = 2;  % Startindex Zeit-Iterationsschleife
  case {'AB','ABC'}
    % Adams-Bashforth (Aufgabe 2) / Adams-Bashforth-Moulton (Aufgabe 3): Anlaufstueck mit zwei Stuetzstellen aus analytischer Loesung
    % t  = ...
    % yh = ...
    % y_numerisch(:,2) = ...    % Anlaufstueck beruecksichtigen
    n_start = 3;  % Startindex Zeit-Iterationsschleife
end

% Zeit-Iterationsschleife .................................................
for n = n_start:N+1
  % Berechnung Loesungsinkrement
  switch integrator
    case 'HN'
       dy = Heun(t,yh,h,m,D,k);
       t = t+h;
       yh = yh+dy;
    case 'AB'
      % dy = AdamsBashforth(...);
    case 'ABM'
      % dy = AdamsBashforthMoulton(...);
  end
    
  % Lauf- und Loesungsvariable inkrementieren


  % Werte in Lösungsvektor speichern
  y_numerisch(:,n) = [t;yh'];
end % .....................................................................

% quadratischen Fehler berechnen und ausgeben
L2Fehler = 0;
for n = n_start:N+1 % verfahrensabhaengigen Startindex beachten 
  % Fehler aufsummieren
  % y_ref = y1_exakt(..);
  % L2Fehler = L2Fehler + ...

  y_ref = y1_exakt_schwingfall(y_numerisch(1,n),y0,m,k,D);
  L2Fehler = L2Fehler + (y_numerisch(2,n)-y_ref)^2;

end

L2Fehler = sqrt( L2Fehler / ((N+1)-(n_start-1))) % Semikolon auskommentieren, um Fehler auszugeben

% analytische Referenzloesung
Nplot = 10*N+1;                 % Anzahl Stuetzstellen fuer Plot: hoehere Sampling-Rate als fuer numerische Loesung
y_analytisch = zeros(2,Nplot);  % erste Zeile = Stuetzstellen t, zweite Zeile = Loesungswerte y

% Fuellen des analytischen Loesungsvektors
y_analytisch(1,:) = linspace(t0,t_end,Nplot);
y_analytisch(2,:) = y1_exakt_schwingfall(y_analytisch(1,:),y0,m,k,D);

% analytische und numerische Loesung plotten
figure()
plot(y_numerisch(1,:),y_numerisch(2,:),'rs',y_analytisch(1,:),y_analytisch(2,:),'b')
title(sprintf('N = %D',N))

end % Funktion AWP_Loeser()


%% UNTERFUNKTIONEN
% =========================================================================

% rechte Seite (RHS) der DGL
function [RHS] = f(t,y,m,D,k)
       z1 = y(1);
       z2 = y(2);
       RHS = [z2 (1/m*(-k*z1 -D*z2))];
  end 

% Loesungsinkrement dy der Integrationsverfahren --------------------------
% Heun-Verfahren
function [dy] = Heun(t,yh,h,m,D,k)
     %D
   y_numerisch = yh + h * f(t,yh,m,D,k);
   dy = (h/2) * (f(t,yh,m,D,k) + f(t + h, y_numerisch, m, D ,k));
end

% Adams-Bashforth-Verfahren 2. Ordnung
function [dy] = AdamsBashforth()
    % dy = ...
end

% Adams-Moulton-Verfahren 3. Ordnung mit Adams-Bashforth-Praediktor
function [dy] = AdamsBashforthMoulton()
    % dy = ...
end


% Analytische Loesungen ---------------------------------------------------

% analytische Loesungen y(t) ----------------------------------------------
% Schwingfall
function [output] = y1_exakt_schwingfall(t,y0,m,k,D)
  delta   = D/(2*m);
  omega   = sqrt(k/m-delta^2);
  output  = y0.*exp(-delta.*t).*(cos(omega.*t)+delta./omega*sin(omega.*t));
end

% aperiodischer Grenzfall
function [output] = y1_exakt_aperiodisch(t,y0,m,k,D)
  delta   = D/(2*m);
  output  = (y0+delta.*y0.*t).*exp(-delta.*t);
end

% Kriechfall
function [output] = y1_exakt_kriechfall(t,y0,m,k,D)
  delta   = D/(2*m);
  omega   = sqrt(delta^2-k/m);
  output = y0.*1./(2.*omega).*((omega+delta).*exp(omega.*t)+(omega-delta).*exp(-omega.*t)).*exp(-delta.*t);
end

% analytische Loesungen y(t)' ---------------------------------------------
% Schwingfall
function [output] = y2_exakt_schwingfall(t,y0,m,k,D)
  delta   = D/(2*m);
  omega   = sqrt(k/m-delta^2);
  output  = -y0.*delta.*exp(-delta.*t).*(cos(omega.*t)+delta./omega.*sin(omega.*t)) ...
            +y0.*exp(-delta*t).*(-omega.*sin(omega.*t)+delta*cos(omega.*t));
end

% aperiodischer Grenzfall
function [output] = y2_exakt_aperiodisch(t,y0,m,k,D)
  delta   = D/(2*m);
  output  = (delta.*y0).*exp(-delta.*t)-(y0+delta.*y0.*t).*delta.*exp(-delta.*t);
end

% Kriechfall
function [output] = y2_exakt_kriechfall(t,y0,m,k,D)
  delta   = D/(2*m);
  omega   = sqrt(delta^2-k/m);
  output  = +y0.*1./(2.*omega).*(omega*(omega+delta).*exp(omega.*t)-omega.*(omega-delta).*exp(-omega.*t)).*exp(-delta.*t) ...
            -delta.*y0.*1./(2.*omega).*((omega+delta).*exp(omega.*t)+(omega-delta).*exp(-omega.*t)).*exp(-delta.*t);
end