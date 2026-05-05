# Controle/Ferramentas/DBC.jl
#
# Otimização SOS via Recurrent Dissipativity-Based Conditions (DBCs)
# Referência: Madeira & Machado, MICNON 2024
#
# USO:
#   K = otimizarDBC(f, g, h, vars_x, u_var)
#   K = otimizarDBC(f, g, h, vars_x, u_var; imax=100, dim_u=1)

using JuMP
using SumOfSquares
using DynamicPolynomials
using LinearAlgebra
using Clarabel
# using COSMO
using Printf

"""
        otimizarDBC(f, g, h, vars_x, u_var;
                    imax       = 100,
                    grau_V     = 1:2,
                    grau_T     = 1:4,
                    epsilon_R  = 1e-4,
                    beta       = 1e-8,
                    verbose    = true) -> Matrix

    Busca um controlador estático de saída u = K*h(x) via otimização SOS com
    Recurrent Dissipativity-Based Inequalities (DBIs).

    Argumentos obrigatórios:
        - f       : Vetor de polinômios [n×1] representando a dinâmica f(x)
        - g       : Vetor/matriz de polinômios ou reais [n×m] — matriz de entrada
        - h       : Vetor de polinômios [p×1] — saída fictícia y = h(x)
        - vars_x  : Vetor de variáveis de estado (@polyvar)
        - u_var   : Variável simbólica de entrada (@polyvar) — escalar ou vetor

    Argumentos opcionais:
        - imax      : Número máximo de iterações (padrão: 100)
        - grau_V    : Faixa de graus dos monômios de V (padrão: 1:2)
        - grau_T    : Faixa de graus dos monômios de T (padrão: 1:4)
        - epsilon_R : Margem de positividade em R (padrão: 1e-4)
        - beta      : Margem de positividade de V e T (padrão: 1e-8)
        - verbose   : Imprime progresso (padrão: true)

    Retorna:
    K :: Matrix{Float64} — ganho do controlador (dim_u × dim_y)
    tal que u = K * h(x) estabiliza o sistema em malha fechada.
"""
function otimizarDBC(f, g, h, vars_x, u_var;
                     imax      ::Int     = 100,
                     grau_V             = 1:2,
                     grau_T             = 1:4,
                     epsilon_R ::Float64 = 1e-4,
                     beta      ::Float64 = 1e-8,
                     verbose   ::Bool    = true)

    # ── Inferência de dimensões ──────────────────────────────────────────
    dim_x = length(vars_x)
    dim_y = length(h)
    dim_u = u_var isa AbstractVector ? length(u_var) : 1

    u_vec = u_var isa AbstractVector ? u_var : [u_var]

    # Calcula a dinâmica afim: f(x) + g(x)*u
    # g pode ser Matrix{Float64} ou vetor — normaliza para matriz [n × m]
    G = _normalizar_g(g, dim_x, dim_u)

    verbose && _banner("Otimização SOS (Recurrent DBIs)", dim_x, dim_y, dim_u)

    # ── Monômios base ────────────────────────────────────────────────────
    Z_v = monomials(vars_x, grau_V)
    Z_t = monomials(vars_x, grau_T)

    # ── Iteração 0 ───────────────────────────────────────────────────────
    Q0_val, S0_val, R0_val = _iteracao0(
        f, G, h, vars_x, u_vec, Z_v, Z_t, dim_y, dim_u,
        epsilon_R, beta, verbose
    )

    Q_final, S_final, R_final = Q0_val, S0_val, R0_val

    # ── Loop iterativo ───────────────────────────────────────────────────
    for i in 1:imax
        Q_val, S_val, R_val, ok = _iteracao_k(
            f, G, h, vars_x, u_vec, Z_v, Z_t, dim_y, dim_u,
            Q0_val, S0_val, R0_val, epsilon_R, beta, verbose, i
        )

        if !ok
            break
        end

        Δi     = S_val * inv(R_val) * S_val' - Q_val
        eig_Δi = eigvals(Δi)[1]

        if eig_Δi >= 0
            Q_final, S_final, R_final = Q_val, S_val, R_val
            verbose && println("  ✓ Solução estabilizante encontrada!!")
            break
        end

        Q0_val, S0_val, R0_val = Q_val, S_val, R_val
    end

    # ── Ganho final ──────────────────────────────────────────────────────
    K = -inv(R_final) * S_final'

    if verbose
        println("\n  Ganho do Controlador K:")
        for i in axes(K, 1)
            @printf("    K[%d,:] = %s\n", i, join([@sprintf("%.4f", v) for v in K[i,:]], "  "))
        end
        println("═"^52, "\n")
    end

    return K
end

# =========================================================================
# FUNÇÕES INTERNAS
# =========================================================================

# Monta a dinâmica afim f + G*u como vetor de polinômios
function _dinamica_afim(f, G, u_vec)
    return [f[i] + sum(G[i, j] * u_vec[j] for j in axes(G, 2)) for i in axes(G, 1)]
end

# Calcula ∇V · (f + G*u) — Lie derivative
function _lie(V_var, vars_x, f_afim)
    grad = differentiate(V_var, vars_x)
    return dot(grad, f_afim)
end

# Função de supply rate r(y, u) = y'Qy + 2y'Su + u'Ru
function _supply(y, u_vec, Q_var, S_var, R_var)
    return dot(y, Q_var * y) +
           2 * dot(y, S_var * u_vec) +
           dot(u_vec, R_var * u_vec)
end

# Normaliza g para Matrix{Any} de tamanho [dim_x × dim_u]
function _normalizar_g(g, dim_x, dim_u)
    if g isa AbstractMatrix
        return g
    elseif g isa AbstractVector && dim_u == 1
        # Coluna única → matriz n×1
        return reshape(Any[gi for gi in g], dim_x, 1)
    else
        error("Formato de g não reconhecido. Passe uma Matrix [n×m] ou Vector [n] para dim_u=1.")
    end
end

# Iteração 0 — sem restrições recorrentes
function _iteracao0(f, G, h, vars_x, u_vec, Z_v, Z_t, dim_y, dim_u,
                    epsilon_R, beta, verbose)

    # model = SOSModel(Clarabel.Optimizer)
     model = SOSModel(optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter"                         => 500,
        "time_limit"                       => 180.0,
        "verbose"                          => false,

        # Tolerâncias relaxadas — combate o NUMERICAL_ERROR
        "tol_gap_abs"                      => 1e-6,
        "tol_gap_rel"                      => 1e-6,
        "tol_feas"                         => 1e-6,
        "tol_infeas_abs"                   => 1e-6,
        "tol_infeas_rel"                   => 1e-6,

        # Tolerâncias reduzidas (usadas quando o solver está perto do limite)
        "reduced_tol_gap_abs"              => 1e-4,
        "reduced_tol_gap_rel"              => 1e-4,
        "reduced_tol_feas"                 => 1e-4,

        # Equilibração — reescala interna das matrizes SDP
        "equilibrate_enable"               => true,
        "equilibrate_max_iter"             => 50,     # padrão: 10
        "equilibrate_min_scaling"          => 1e-6,
        "equilibrate_max_scaling"          => 1e6,

        # Regularização estática — estabiliza matrizes quase singulares
        "static_regularization_enable"     => true,
        "static_regularization_constant"   => 1e-7,  # padrão: 1e-8, sobe levemente

        # Regularização dinâmica — ativa quando autovalores ficam muito pequenos
        "dynamic_regularization_enable"    => true,
        "dynamic_regularization_eps"       => 1e-12,
        "dynamic_regularization_delta"     => 1e-6,  # padrão: 2e-7, sobe levemente

        # Refinamento iterativo — melhora a solução após cada passo
        "iterative_refinement_enable"      => true,
        "iterative_refinement_max_iter"    => 20,    # padrão: 10
        "iterative_refinement_reltol"      => 1e-12,
        "iterative_refinement_abstol"      => 1e-11,
    ))
    set_silent(model)

    @variable(model, V0, SOSPoly(Z_v))
    @variable(model, T0, SOSPoly(Z_t))
    @variable(model, Q0[1:dim_y, 1:dim_y], Symmetric)
    @variable(model, S0[1:dim_y, 1:dim_u])
    @variable(model, R0[1:dim_u, 1:dim_u], Symmetric)

    # Positividade de R
    for j in 1:dim_u
        @constraint(model, R0[j, j] >= epsilon_R)
    end

    # Positividade de V e T
    x_sq = sum(v^2 for v in vars_x)
    @constraint(model, V0 - beta * x_sq        in SOSCone())
    @constraint(model, T0 - beta * x_sq^2      in SOSCone())

    # Dissipatividade
    f_afim = _dinamica_afim(f, G, u_vec)
    Vdot   = _lie(V0, vars_x, f_afim)
    r_yu   = _supply(h, u_vec, Q0, S0, R0)
    @constraint(model, -Vdot - T0 + r_yu in SOSCone())

    optimize!(model)

    if termination_status(model) ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        error("Solver falhou na iteração 0. Tente ajustar grau_V, grau_T ou epsilon_R.")
    end

    Q0_val = value.(Q0)
    S0_val = value.(S0)
    R0_val = value.(R0)

    if verbose
        Δ0 = S0_val * inv(R0_val) * S0_val' - Q0_val
        @printf("  Iteração   0 → EIG(Δ) = %+.5f\n", eigvals(Δ0)[1])
    end

    return Q0_val, S0_val, R0_val
end

# Iteração k — com restrições recorrentes (LMIs)
function _iteracao_k(f, G, h, vars_x, u_vec, Z_v, Z_t, dim_y, dim_u,
                     Q0_val, S0_val, R0_val, epsilon_R, beta, verbose, i)

    # model = SOSModel(Clarabel.Optimizer)
    model = SOSModel(optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter"                         => 500,
        "time_limit"                       => 180.0,
        "verbose"                          => false,

        # Tolerâncias relaxadas — combate o NUMERICAL_ERROR
        "tol_gap_abs"                      => 1e-6,
        "tol_gap_rel"                      => 1e-6,
        "tol_feas"                         => 1e-6,
        "tol_infeas_abs"                   => 1e-6,
        "tol_infeas_rel"                   => 1e-6,

        # Tolerâncias reduzidas (usadas quando o solver está perto do limite)
        "reduced_tol_gap_abs"              => 1e-4,
        "reduced_tol_gap_rel"              => 1e-4,
        "reduced_tol_feas"                 => 1e-4,

        # Equilibração — reescala interna das matrizes SDP
        "equilibrate_enable"               => true,
        "equilibrate_max_iter"             => 50,     # padrão: 10
        "equilibrate_min_scaling"          => 1e-6,
        "equilibrate_max_scaling"          => 1e6,

        # Regularização estática — estabiliza matrizes quase singulares
        "static_regularization_enable"     => true,
        "static_regularization_constant"   => 1e-7,  # padrão: 1e-8, sobe levemente

        # Regularização dinâmica — ativa quando autovalores ficam muito pequenos
        "dynamic_regularization_enable"    => true,
        "dynamic_regularization_eps"       => 1e-12,
        "dynamic_regularization_delta"     => 1e-6,  # padrão: 2e-7, sobe levemente

        # Refinamento iterativo — melhora a solução após cada passo
        "iterative_refinement_enable"      => true,
        "iterative_refinement_max_iter"    => 20,    # padrão: 10
        "iterative_refinement_reltol"      => 1e-12,
        "iterative_refinement_abstol"      => 1e-11,
    ))
    set_silent(model)

    @variable(model, V, SOSPoly(Z_v))
    @variable(model, T, SOSPoly(Z_t))
    @variable(model, Q[1:dim_y, 1:dim_y], Symmetric)
    @variable(model, S[1:dim_y, 1:dim_u])
    @variable(model, R[1:dim_u, 1:dim_u], Symmetric)

    for j in 1:dim_u
        @constraint(model, R[j, j] >= epsilon_R)
    end

    x_sq = sum(v^2 for v in vars_x)
    @constraint(model, V - beta * x_sq   in SOSCone())
    @constraint(model, T - beta * x_sq^2 in SOSCone())

    f_afim = _dinamica_afim(f, G, u_vec)
    Vdot   = _lie(V, vars_x, f_afim)
    r_yu   = _supply(h, u_vec, Q, S, R)
    @constraint(model, -Vdot - T + r_yu in SOSCone())

    # Restrições recorrentes
    M = S*inv(R0_val)*S0_val' +
        S0_val*inv(R0_val)*S' -
        2*S0_val*inv(R0_val)*S0_val' +
        Q0_val - Q

    @constraint(model, Symmetric(M)          in PSDCone())
    @constraint(model, Symmetric(R0_val - R) in PSDCone())

    optimize!(model)

    status = termination_status(model)
    if status ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        verbose && @printf("  [!] Solver falhou na iteração %d (status: %s). Parando.\n", i, status)
        return Q0_val, S0_val, R0_val, false
    end

    Q_val = value.(Q)
    S_val = value.(S)
    R_val = value.(R)

    Δi     = S_val * inv(R_val) * S_val' - Q_val
    eig_Δi = eigvals(Δi)[1]
    verbose && @printf("  Iteração %3d → EIG(Δ) = %+.5f\n", i, eig_Δi)

    return Q_val, S_val, R_val, true
end

# Banner de abertura
function _banner(titulo, dim_x, dim_y, dim_u)
    println("═"^52)
    println("  ", titulo)
    @printf("  Estados: %d  |  Saída: %d  |  Entradas: %d\n", dim_x, dim_y, dim_u)
    println("═"^52)
end