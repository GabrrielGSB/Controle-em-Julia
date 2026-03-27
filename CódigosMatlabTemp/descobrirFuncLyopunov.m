%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOSDEMO2 - Lyapunov Function Search% Section 3.2 of SOSTOOLS User’s Manual

% Exemplo baseado em SOSDEMO2

clc, clear all

syms x1 x2 x3;
vars = [x1; x2; x3];

% u(x) - State feedback control law
u = -3.6345*x1^3 + 4.4439*x1^2*x2 - 7.5113*x1*x2^2 - 3.5452*x2^3; 

% Constructing the vector field dx/dt = f + gu, g= [0; 1]
f = [2*x1^3+x1^2*x2-6*x1*x2^2+5*x2^3; 0] + [0; 1]*u;

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);

% =============================================
% The Lyapunov function V(x):
[prog,V] = sospolyvar(prog,[monomials([x1,x2],[2,3,4,5,6])]);

% =============================================
% Next, define SOSP constraints

% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
prog = sosineq(prog,V-10^-2*(x1^6+x2^6));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2));
prog = sosineq(prog,expr);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [Q1,Z1] = findsos(SOLV)
%  simplify(  SOLV - simplify(Z1'*Q1*Z1)  )