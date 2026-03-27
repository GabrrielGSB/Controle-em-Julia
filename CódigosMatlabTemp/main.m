clc, clear all

syms x1 x2 x3;
vars = [x1; x2];

f = [2*x1^3+x1^2*x2-6*x1*x2^2+5*x2^3; 0];
g = [0; 1];

V0 = 1*vars'*vars;

[u, lamb] = calcularControleSOS(f, g, V0, vars);