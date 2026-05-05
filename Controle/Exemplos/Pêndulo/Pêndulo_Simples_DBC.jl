# Controle/Exemplos/Pêndulo/Pêndulo_DBC.jl
#
# Estabilização do Pêndulo Simples via
# Dissipativity-Based Conditions (DBC) com otimização SOS.
#
# Modelo polinomial obtido via aproximação de Taylor:
#   sin(θ) ≈ θ - θ³/6
#
# Adimensionalização: τ = ω₀·t, onde ω₀ = √(g/L)
#
# Lei de controle sintetizada (no sistema adimensional):
#   u_nd = K_nd · h(x),  h(x) = [x₁, x₂]
#   u_dim = u_nd · (m·L²·ω₀)   ← torque físico aplicado

using DynamicPolynomials

using Printf
# =========================================================
# INCLUDES
    include("../../SistemasLTI/Pêndulo.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
    include("../../Ferramentas/DBC.jl")
# =========================================================

# =========================================================
# 1. PARÂMETROS FÍSICOS
    L_val = 1.0
    m_val = 1.0
    b_val = 0.1
    g_val = 9.81

    I_val = m_val * L_val^2

    x0    = [π/4, 0.0]
    tspan = (0.0, 20.0)
# =========================================================

# =========================================================
# 2. ADIMENSIONALIZAÇÃO
    # Reescala de tempo: τ = ω₀·t  →  d/dt = ω₀ · d/dτ
    #
    # Sistema adimensional resultante:
    #   dx₁/dτ =  x₂
    #   dx₂/dτ = -(x₁ - x₁³/6) - β_nd·x₂  +  u_nd
    #
    # onde:
    #   β_nd  = b / (m·L²·ω₀)    — amortecimento adimensional
    #   u_nd  = u_dim / (m·L²·ω₀²)  — entrada adimensional
    #
    # Ao recuperar o ganho físico:
    #   u_dim = u_nd · (m·L²·ω₀²)
    #   K_dim = K_nd / (m·L²·ω₀²)   ← divide porque u_dim = K_dim·h(x)

    ω0   = sqrt(g_val / L_val)                # frequência natural ≈ 3.13 rad/s
    β_nd = b_val / (m_val * L_val^2 * ω0)    # ≈ 0.032  — O(1) ✓

    @printf("Adimensionalização:\n")
    @printf("  ω₀   = %.4f rad/s\n", ω0)
    @printf("  β_nd = %.4f  (era b/I = %.4f)\n\n", β_nd, b_val/I_val)
# =========================================================

# =========================================================
# 3. MODELO POLINOMIAL ADIMENSIONAL
    @polyvar x1 x2 u_sym

    f = [ x2, -(x1 - x1^3 / 6) - β_nd * x2 ]
    g = [0.0, 1.0]  
    h = [x1, x2]
# =========================================================

# =========================================================
# 4. SÍNTESE DO CONTROLADOR VIA DBC (no espaço adimensional)
    println("Sintetizando controlador DBC para o Pêndulo Simples...\n")

    K_nd = otimizarDBC(
        f, g, h,
        [x1, x2], u_sym;
        grau_V    = 1:4,
        grau_T    = 1:4,
        epsilon_R = 1e-3,
        beta      = 1e-6,
        imax      = 1000,
        verbose   = true
    )
# =========================================================

# =========================================================
# 5. DESFAZENDO A ESCALA DO GANHO
    # u_dim = u_nd · (m·L²·ω₀²)
    # u_nd  = K_nd · [x₁, x₂]
    # u_dim = K_nd · (m·L²·ω₀²) · [x₁, x₂]
    #       = K_dim · [x₁, x₂]

    escala_u = m_val * L_val^2 * ω0^2    # = m·L²·ω₀² = m·g·L  (= I·ω₀²)

    K1 = K_nd[1, 1] * escala_u
    K2 = K_nd[1, 2] * escala_u

    @printf("\nGanhos no espaço adimensional : K_nd = [%.4f  %.4f]\n", K_nd[1,1], K_nd[1,2])
    @printf("Ganhos físicos (torque [N·m]) : K₁   = %.4f   K₂ = %.4f\n", K1, K2)
    @printf("Lei de controle: u(x) = %.4f·θ + %.4f·ω\n\n", K1, K2)
# =========================================================

# =========================================================
# 6. STRUCT E DINÂMICAS PARA SIMULAÇÃO
    struct PenduloDBCParams
        L               ::Float64
        m               ::Float64
        b               ::Float64
        K1              ::Float64
        K2              ::Float64
        estadosIniciais ::Vector{Float64}
        numEstados      ::Int
    end

    function PenduloDBCParams(; L, m, b, K1, K2, estadosIniciais)
        PenduloDBCParams(L, m, b, K1, K2, estadosIniciais, 2)
    end

    # Dinâmica exata em malha fechada — sin(x₁) sem aproximação
    function pendulo_dbc_mf!(dx, x, p::PenduloDBCParams, t)
        I = p.m * p.L^2
        u = p.K1 * x[1] + p.K2 * x[2]

        dx[1] = x[2]
        dx[2] = (-p.m * 9.81 * p.L * sin(x[1]) - p.b * x[2] + u) / I
    end

    # Dinâmica em malha aberta
    function pendulo_aberta!(dx, x, p::PenduloDBCParams, t)
        I = p.m * p.L^2
        dx[1] = x[2]
        dx[2] = (-p.m * 9.81 * p.L * sin(x[1]) - p.b * x[2]) / I
    end
# =========================================================

# =========================================================
# 7. SIMULAÇÃO
    params = PenduloDBCParams(
        L  = L_val,
        m  = m_val,
        b  = b_val,
        K1 = K1,
        K2 = K2,
        estadosIniciais = x0
    )

    sol_aberta  = resolverSistema(pendulo_aberta!,   x0, tspan, params; resolucao=0.005)
    sol_fechada = resolverSistema(pendulo_dbc_mf!,   x0, tspan, params; resolucao=0.005)
# =========================================================

# =========================================================
# 8. VISUALIZAÇÃO
    plotarNoTempo(
        [sol_aberta, sol_fechada];
        nomes   = ["Malha Aberta", "DBC (SOS)"],
        titulo  = "Pêndulo — Comparação Temporal: Sem Controle vs. DBC",
        estados = 1:2
    )

    plotarRetratoFase(
        [sol_aberta, sol_fechada];
        estados = (1, 2),
        nomes   = ["Malha Aberta", "DBC (SOS)"],
        titulo  = "Pêndulo — Retrato de Fase (θ × ω)"
    )

    plotarRetratoFaseCompleto(
        pendulo_dbc_mf!, params;
        limiteEixo         = 2.0,
        densidadeSetas     = 18,
        raiosIniciais      = [0.3, 0.8, 1.4],
        trajetoriasPorAnel = 3,
        tempoMaximo        = 15.0,
        titulo             = "Pêndulo DBC — Campo Vetorial em Malha Fechada"
    )
# =========================================================

# =========================================================
# 9. ANÁLISE DE PERFORMANCE
    m_theta = analisarPerformance(sol_fechada;
        referencia = 0.0,
        idx_estado = 1,
        banda      = 0.05
    )

    m_omega = analisarPerformance(sol_fechada;
        referencia = 0.0,
        idx_estado = 2,
        banda      = 0.05
    )

    imprimirRelatorio(m_theta, nome = "Pêndulo DBC — Ângulo θ")
    imprimirRelatorio(m_omega, nome = "Pêndulo DBC — Velocidade angular ω")
# =========================================================

# =========================================================
# 10. COMPARAÇÃO COM PID
    include("../../Controladores/PID.jl")
    include("../../Abstrações/MalhaFechada.jl")

    pendulo_planta = Pendulo(
        comprimento     = L_val,
        massa           = m_val,
        coefAtrito      = b_val,
        estadosIniciais = x0
    )

    pid_ref = PID(Kp = 30.0, Ki = 5.0, Kd = 8.0)
    sys_pid = conectar(pendulo_planta, pid_ref, 0.0)
    x0_pid  = condicoesIniciais(sys_pid, x0)
    sol_pid = resolverSistema(sys_pid, x0_pid, tspan; resolucao=0.005)

    println("\n─── Comparação de Ganhos ─────────────────────────────────")
    @printf("  K DBC  →  K₁ = %+.4f   K₂ = %+.4f\n", K1, K2)
    @printf("  PID    →  Kp = %.1f   Ki = %.1f   Kd = %.1f\n",
            pid_ref.Kp, pid_ref.Ki, pid_ref.Kd)
    println("──────────────────────────────────────────────────────────\n")

    compararControladores(
        [sol_pid, sol_fechada],
        ["PID (referência)", "DBC (SOS)"];
        referencia = 0.0,
        idx_estado = 1,
        banda      = 0.05
    )
# =========================================================