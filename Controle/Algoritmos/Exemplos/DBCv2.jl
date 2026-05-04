# Controle/Exemplos/DBC_Dissipativo.jl
#
# Implementação SOS - Recurrent DBIs (Dissipativity-Based Conditions)
# Referência: Madeira e Machado, MICNON 2024
#             Valmorbida and Papachristodoulou (Exemplo 2)
#

# =========================================================
# INCLUDES
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
# =========================================================

using JuMP
using SumOfSquares
using DynamicPolynomials
using LinearAlgebra
using Clarabel
using OrdinaryDiffEq
using DiffEqCallbacks
using Printf

# =========================================================================
# 1. DEFINIÇÃO DA PLANTA POLINOMIAL
# =========================================================================

#=
    Sistema polinomial do Exemplo 2 (Valmorbida & Papachristodoulou):
        ẋ₁ = 2x₁³ + x₁²x₂ - 6x₁x₂² + 5x₂³
        ẋ₂ = u
    
    Saída fictícia para o projeto de dissipatividade:
        h(x) = x₁³ + x₁²x₂ + x₁x₂² + x₂³
    
    Lei de controle (malha fechada):
        u = K * h(x)
=#

function sistema_dbc!(dx, x, p, t; u=0.0)
    x1, x2 = x[1], x[2]

    dx[1] = 2*x1^3 + (x1^2)*x2 - 6*x1*x2^2 + 5*x2^3
    dx[2] = u
end

# Wrapper com ganho K embutido — compatível com resolverSistema
struct SistemaDBCParams
    K               ::Float64
    estadosIniciais ::Vector{Float64}
    numEstados      ::Int
end

function SistemaDBCParams(; K, estadosIniciais)
    SistemaDBCParams(K, estadosIniciais, 2)
end

function sistema_dbc_mf!(dx, x, p::SistemaDBCParams, t)
    x1, x2 = x[1], x[2]
    h = x1^3 + x1^2*x2 + x1*x2^2 + x2^3
    u = p.K * h
    sistema_dbc!(dx, x, p, t; u=u)
end

# =========================================================================
# 2. OTIMIZAÇÃO SOS — RECURRENT DBIs
# =========================================================================

function otimizarDBC(; imax=100, epsilon_R=1e-4, beta=1e-8, verbose=true)
    @polyvar x1 x2 u
    x_vars = [x1, x2]

    f = [2*x1^3 + (x1^2)*x2 - 6*x1*x2^2 + 5*(x2^3), 0]
    g = [0, 1]
    h_poly = x1^3 + x1^2*x2 + x1*x2^2 + x2^3
    y = [h_poly]

    # ── Iteração 0 ──────────────────────────────────────────────────────
    verbose && println("═"^50)
    verbose && println("  Iniciando Otimização SOS (Recurrent DBIs)")
    verbose && println("═"^50)

    Z_v = monomials(x_vars, 1:2)

    model0 = SOSModel(Clarabel.Optimizer)
    set_silent(model0)

    @variable(model0, V0, SOSPoly(Z_v))
    @variable(model0, T0, SOSPoly(Z_v))
    @variable(model0, Q0[1:1, 1:1], Symmetric)
    @variable(model0, S0[1:1, 1:1])
    @variable(model0, R0[1:1, 1:1], Symmetric)

    @constraint(model0, R0[1,1] >= epsilon_R)
    @constraint(model0, V0 - beta*(x1^2 + x2^2)   in SOSCone())
    @constraint(model0, T0 - beta*(x1^2 + x2^2)^2 in SOSCone())

    gradV0 = differentiate(V0, x_vars)
    Vdot0  = dot(gradV0, f + g*u)
    r_yu0  = dot(y, Q0*y) + 2*dot(y, S0*[u]) + dot([u], R0*[u])
    @constraint(model0, -Vdot0 - T0 + r_yu0 in SOSCone())

    optimize!(model0)

    if termination_status(model0) ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        error("Solver falhou na iteração 0.")
    end

    Q0_val = value.(Q0)
    S0_val = value.(S0)
    R0_val = value.(R0)

    Δ0 = S0_val*inv(R0_val)*S0_val' - Q0_val
    verbose && @printf("  Iteração 0 → EIG(Δ) = %.5f\n", eigvals(Δ0)[1])

    Q_final, S_final, R_final = Q0_val, S0_val, R0_val

    # ── Loop iterativo ───────────────────────────────────────────────────
    for i in 1:imax
        model = SOSModel(Clarabel.Optimizer)
        set_silent(model)

        @variable(model, V, SOSPoly(Z_v))
        @variable(model, T, SOSPoly(Z_v))
        @variable(model, Q[1:1, 1:1], Symmetric)
        @variable(model, S[1:1, 1:1])
        @variable(model, R[1:1, 1:1], Symmetric)

        @constraint(model, R[1,1] >= epsilon_R)
        @constraint(model, V - beta*(x1^2 + x2^2)   in SOSCone())
        @constraint(model, T - beta*(x1^2 + x2^2)^2 in SOSCone())

        gradV = differentiate(V, x_vars)
        Vdot  = dot(gradV, f + g*u)
        r_yu  = dot(y, Q*y) + 2*dot(y, S*[u]) + dot([u], R*[u])
        @constraint(model, -Vdot - T + r_yu in SOSCone())

        # Restrições recorrentes (LMIs)
        M = S*inv(R0_val)*S0_val' + S0_val*inv(R0_val)*S' - 2*S0_val*inv(R0_val)*S0_val' + Q0_val - Q
        @constraint(model, Symmetric(M)          in PSDCone())
        @constraint(model, Symmetric(R0_val - R) in PSDCone())

        optimize!(model)

        if termination_status(model) ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
            verbose && println("  [!] Solver falhou na iteração $i. Interrompendo.")
            break
        end

        Q_val = value.(Q)
        S_val = value.(S)
        R_val = value.(R)

        Δi     = S_val*inv(R_val)*S_val' - Q_val
        eig_Δi = eigvals(Δi)[1]
        verbose && @printf("  Iteração %3d → EIG(Δ) = %.5f\n", i, eig_Δi)

        if eig_Δi >= 0
            Q_final, S_final, R_final = Q_val, S_val, R_val
            verbose && println("  ✓ Solução estabilizante encontrada!")
            break
        end

        Q0_val, S0_val, R0_val = Q_val, S_val, R_val
    end

    K = -inv(R_final) * S_final'
    verbose && @printf("\n  Ganho do Controlador K = %.4f\n", K[1,1])
    verbose && println("═"^50, "\n")

    return K[1,1]
end

# =========================================================================
# 3. EXECUÇÃO
# =========================================================================

using Printf

K_calculado = otimizarDBC(imax=100, verbose=true)

x0       = [0.5, -0.5]
tspan    = (0.0, 15.0)

# Parâmetros para cada cenário
p_fechada = SistemaDBCParams(K=K_calculado,  estadosIniciais=x0)
p_aberta  = SistemaDBCParams(K=0.0,          estadosIniciais=x0)

# ── Malha fechada ────────────────────────────────────────────────────────
sol_fechada = resolverSistema(sistema_dbc_mf!, x0, tspan, p_fechada,
                              resolucao=0.005)

# ── Malha aberta (com callback de segurança) ─────────────────────────────
#   O solver é chamado diretamente aqui pois precisamos de um callback
#   personalizado que o resolverSistema padrão não expõe.
prob_aberta = ODEProblem(sistema_dbc_mf!, x0, tspan, p_aberta)
cb_explosao = DiscreteCallback(
    (u, t, i) -> norm(u) > 100.0,
    terminate!
)
sol_aberta = solve(prob_aberta, Tsit5(),
                   callback=cb_explosao,
                   reltol=1e-6, abstol=1e-6,
                   saveat=0.005)

# =========================================================================
# 4. VISUALIZAÇÃO — PADRÃO DO FRAMEWORK
# =========================================================================

# ── Estados ao longo do tempo (apenas malha fechada é estável) ───────────
plotarNoTempo(sol_fechada,
              titulo  = "DBC — Estados em Malha Fechada (K = $(round(K_calculado, digits=4)))",
              estados = 1:2)

# ── Retrato de fase comparativo: MF vs MA ───────────────────────────────
plotarRetratoFase([sol_aberta, sol_fechada],
                  estados = (1, 2),
                  nomes   = ["Malha Aberta (diverge)", "Malha Fechada (DBC)"],
                  titulo  = "DBC — Retrato de Fase: Malha Aberta vs Fechada")

                  
function sistema_dbc_campo!(dx, x, p, t)
    x1, x2 = x[1], x[2]
    h = x1^3 + x1^2*x2 + x1*x2^2 + x2^3
    u = p.K * h
    dx[1] = 2*x1^3 + (x1^2)*x2 - 6*x1*x2^2 + 5*x2^3
    dx[2] = u
end

plotarRetratoFaseCompleto(sistema_dbc_campo!, p_fechada,
                          limiteEixo         = 3,
                          densidadeSetas     = 18,
                          raiosIniciais      = [0.3, 0.7, 1.2],
                          trajetoriasPorAnel = 8,
                          tempoMaximo        = 5.0,
                          titulo             = "DBC — Campo Vetorial em Malha Fechada")

plotarRetratoFaseCompleto(sistema_dbc!, p_aberta,
                          limiteEixo         = 1,
                          densidadeSetas     = 18,
                          raiosIniciais      = [0.3, 0.7],
                          trajetoriasPorAnel = 3,
                          tempoMaximo        = 1.0,
                          titulo             = "DBC — Campo Vetorial em Malha Fechada")