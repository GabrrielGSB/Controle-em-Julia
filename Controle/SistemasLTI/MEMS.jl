#=
1. SISTEMA MICROELETROMECÂNICO (MEMS)-----------------------------------------------------------------
    Obs: Sistema extraido de Necessary and Sufficient Dissipativity-Based Conditions for Feedback Stabilization
         por Diego de S. Madeira.

    1.1. DESCRIÇÃO:
        Representa o modelo normalizado de um Sistema Microeletromecânico (MEMS) 
        válido localmente em torno da origem. O sistema modela a dinâmica de deflexão 
        mecânica de uma placa acoplada às forças eletrostáticas provenientes de 
        uma carga elétrica.

    1.2. MODELAGEM EM ESPAÇO DE ESTADOS:
        O sistema é não linear e decomposto em três variáveis de estado 
        acopladas:

        1.2.1. Vetor de Variáveis de Estado[cite: 392]: 
            - x[1] : Deflexão da placa
            - x[2] : Velocidade de deflexão
            - x[3] : Carga elétrica

        1.2.2. Vetor de Derivadas (dx) em malha fechada (u = K * x[3]):
            - dx[1] = x2                                                          -> Dinâmica da deflexão.
            - dx[2] = -x1 - 2*ξ*x2 + (2*q̄/3)*x3 + (1/3)*x3^2                      -> Dinâmica da velocidade (aceleração).
            - dx[3] = (q̄/r)*x1 + (1/r)*x1*x3 + ((x_bar - 1.0)/r)*x3 + (2/(3*r))*u -> Dinâmica eletrostática.                               

    1.3. PARÂMETROS[cite: 392]:
        - q_bar (q̄):   Parâmetro nominal de carga
        - xi (ξ):      Fator de amortecimento relativo
        - r_param (r): Constante de acoplamento do sistema
        - x_bar (x̄):   Parâmetro espacial
=#

struct MEMSParams
    q_bar           ::Float64
    xi              ::Float64
    r_param         ::Float64
    x_bar           ::Float64
    K               ::Float64 

    function MEMSParams(; q_bar, xi, r_param, x_bar, K=0.0)
        new(q_bar, xi, r_param, x_bar, K)
    end
end

function MEMS!(dx, x, p, t)
    q̄ = p.q_bar
    ξ = p.xi
    r = p.r_param
    x̄ = p.x_bar
    K = p.K 

    x1, x2, x3 = x[1], x[2], x[3]
    u = K * x3

    dx[1] = x2
    dx[2] = -x1 - 2*ξ*x2 + (2*q̄/3)*x3 + (1/3) * x3^2
    dx[3] = (q̄/r)*x1 + (1/r)*x1*x3 + ((x̄ - 1.0)/r)*x3 + (2/(3*r))*u
end