%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOSDEMO2 - Lyapunov Function Search% Section 3.2 of SOSTOOLS User’s Manual
% Exemplo baseado em SOSDEMO2
clc, clear all
syms x1 x2 x3;
vars = [x1; x2];

V0 = (x1^2 + x2^2);
grad_V0 = [diff(V0,x1) ; diff(V0,x2)];

f = [2*x1^3+x1^2*x2-6*x1*x2^2+5*x2^3; 0];
g = [0; 1];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);
% =============================================

% The Lyapunov function V(x):
% [prog, V] = sospolyvar(prog,[monomials([x1,x2], [1:10])]);
[prog, u]    = sospolyvar(prog,[monomials([x1,x2], [3])]);
[prog, lamb] = sospolyvar(prog,[monomials([x1,x2], [2:6])]);
% =============================================

% Next, define SOSP constraints
% prog = sosineq(prog,V -10^-2*(x1^6+x2^6));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -grad_V0' * (f + g*u) + lamb;
prog = sosineq(prog, expr);

% =============================================
% And call solver
solver_opt.solver = 'mosek';
[prog,info] = sossolve(prog, solver_opt);
% =============================================

% Finally, get solution
SOLu    = sosgetsol(prog, u)
SOLlamb = sosgetsol(prog, lamb)
% SOLV = sosgetsol(prog, V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Q1,Z1] = findsos(SOLu)
% simplify(  SOLV - simplify(Z1'*Q1*Z1)  )
