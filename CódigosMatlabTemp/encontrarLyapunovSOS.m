function SOLV = encontrarLyapunovSOS(f, g, u, vars)
% ENCONTRAR_LYAPUNOV_SOS Busca uma função de Lyapunov V(x) dado um controle u
% Entradas:
%   f    - Vetor do campo vetorial do sistema em malha aberta
%   g    - Matriz/vetor de entrada do sistema
%   u    - Lei de controle polinomial fixa
%   vars - Vetor simbólico das variáveis de estado (ex: [x1; x2; x3])
% Saída:
%   SOLV - Função de Lyapunov polinomial encontrada

    % Extrai as variáveis x1 e x2 do vetor de estados
    % (Baseado no seu script, V depende apenas de x1 e x2)
    x1 = vars(1);
    x2 = vars(2);

    % Constrói o sistema em malha fechada (dx/dt = f + g*u)
    f_cl = f + g * u;

    % Inicializa o programa SOS
    prog = sosprogram(vars);

    % Define a variável polinomial V(x) com monômios de grau 2 a 6
    [prog, V] = sospolyvar(prog, monomials([x1, x2], 2:6));

    % =============================================
    % Restrições SOSP
    % =============================================
    % Restrição 1: V(x) deve ser localmente definida positiva
    prog = sosineq(prog, V - 10^-2*(x1^6 + x2^6));

    % Restrição 2: Derivada de V ao longo do tempo deve ser semi-negativa
    expr = -(diff(V, x1) * f_cl(1) + diff(V, x2) * f_cl(2));
    prog = sosineq(prog, expr);

    % =============================================
    % Configura e chama o solver
    % =============================================
    solver_opt.solver = 'mosek'; % Garantindo o uso do Mosek
    [prog, info] = sossolve(prog, solver_opt);

    % =============================================
    % Retorna a solução
    % =============================================
    if info.pinf == 0 && info.dinf == 0
        SOLV = sosgetsol(prog, V);
        fprintf('Função de Lyapunov V(x) encontrada com sucesso!\n');
        
        
        [Q1, Z1] = findsos(SOLV);
        if ~isempty(Q1)
            disp('Decomposição Z^T * Q * Z encontrada para V.');
        end
    else
        warning('O solver não conseguiu encontrar V (Problema inviável ou erro numérico).');
        SOLV = [];
    end
end