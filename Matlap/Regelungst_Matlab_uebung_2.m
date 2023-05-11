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


