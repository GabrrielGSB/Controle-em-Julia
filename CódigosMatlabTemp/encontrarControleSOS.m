function [SOLu, SOLlamb] = encontrarControleSOS(f, g, V, vars)
% calcular_controle_sos Busca uma lei de controle u(x) e um multiplicador lamb(x)
% Entradas:
%   f    - Vetor simbólico representando a dinâmica do sistema
%   V   - Função de Lyapunov candidata (expressão simbólica)
%   vars - Vetor simbólico das variáveis de estado (ex: [x1; x2])
% Saídas:
%   SOLu    - Polinômio de controle encontrado
%   SOLlamb - Polinômio do multiplicador encontrado

    grad_V = jacobian(V, vars)'; 

    % Inicializa o programa SOS
    prog = sosprogram(vars);

    % Define as variáveis polinomiais u (grau 3) e lamb (graus de 2 a 6)
    [prog, u]    = sospolyvar(prog, monomials(vars, 3));
    [prog, lamb] = sospolyvar(prog, monomials(vars, 2:6));

    % Define a restrição principal: -dV/dx * (f + g*u) + lamb >= 0
    expr = -grad_V'*(f + g*u) + lamb;
    prog = sosineq(prog, expr);

    % Configura e chama o solver 
    solver_opt.solver = 'mosek';
    [prog, info] = sossolve(prog, solver_opt);

    % Extrai os resultados caso a otimização tenha sido bem-sucedida
    if info.pinf == 0 && info.dinf == 0
        SOLu    = sosgetsol(prog, u)
        SOLlamb = sosgetsol(prog, lamb)
        fprintf('Solução encontrada com sucesso!\n');
    else
        warning('O solver não conseguiu encontrar uma solução (Problema inviável ou erro numerico).');
        SOLu    = [];
        SOLlamb = [];
    end
end