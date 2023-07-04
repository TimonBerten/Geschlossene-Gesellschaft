% verwende 1000 Punkte zum Plotten der Funktion
omega = 10;
n_plot = 1000;
x1 = linspace(0,2*pi,n_plot);
y1 = sin(omega*(x1+pi/4));

% n_sample HIER ANPASSEN
% -----------------------
n_sample = 13;
% -----------------------
x2 = linspace(0,2*pi,n_sample);
y2 = sin(omega*(x2+pi/4));

% plotte sin(10*(x+pi/4)) mit verschiedener Anzahl Stuetzstellen
plot(x1,y1,'k',x2,y2,'bs-');