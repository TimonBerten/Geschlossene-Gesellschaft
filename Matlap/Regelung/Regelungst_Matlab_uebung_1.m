clear 
close all 
clc

w = 1;
zeta = 0.25;
w0 = 1;

num = w0^2;
den = [1, 2*zeta*w0, w0^2];

sys_G = tf(num, den);

t = linspace(0, 40, 1000);

u = sin(w *t);

lsim(sys_G, u, t);

grid;

%Hnadische berechnung für Winkel und verschiebung
Ampl = 2;
del_t = -(29.8 - 28.3);

phi = w * del_t;  

A = freqresp(sys_G, j*w);
B = abs(A);
C = angle(A);
C_deg = C * 180/pi;
y_stat = B * sin(w*t + C);
hold on 
plot(t, y_stat, 'r')

%Bode-Diagramm

figure 
bode(sys_G)
grid


%Nyquist-Diagramm

figure
nyquist(sys_G)
axis equal
hold on 
plot(A, 'r*')
grid
