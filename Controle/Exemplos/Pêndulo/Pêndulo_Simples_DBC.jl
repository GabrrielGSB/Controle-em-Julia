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
#   u_ndim = K_ndim · h(x),  h(x) = [x₁, x₂]
#   u_dim = u_ndim · (m·L²·ω₀)   ← torque físico aplicado

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
    L = 1.0
    m = 1.0
    β = 0.1
    g = 9.81

    I = m * L^2

    x0    = [π/4, 0.0]
    tspan = (0.0, 20.0)
# =========================================================

# =========================================================
# 2. ADIMENSIONALIZAÇÃO
    # Reescala de tempo: τ = ω₀·t  →  d/dt = ω₀ · d/dτ
    #
    # Sistema adimensional resultante:
    #   dx₁/dτ =  x₂
    #   dx₂/dτ = -(x₁ - x₁³/6) - β_nd·x₂  +  u_ndim
    #
    # onde:
    #   β_ndim  = β / (m·L²·ω₀)       — amortecimento adimensional
    #   u_ndim  = u_dim / (m·L²·ω₀²)  — entrada adimensional
    #
    # Ao recuperar o ganho físico:
    #   u_dim = u_ndim · (m·L²·ω₀²)
    #   K_dim = K_ndim / (m·L²·ω₀²)   ← divide porque u_dim = K_dim·h(x)

    ω0     = sqrt(g / L)        # frequência natural 
    β_ndim = β / (m * L^2 * ω0)    
# =========================================================

# =========================================================
# 3. MODELO POLINOMIAL ADIMENSIONAL
    @polyvar x1 x2 u_sym

    f = [ x2, -(x1 - x1^3 / 6) - β_ndim * x2 ]
    g = [0.0, 1.0]  
    h = [x1, x2]
# =========================================================

# =========================================================
# 3.1 MODELO POLINOMIAL ADIMENSIONAL (com referência)
    θ_ref = 0

    @polyvar x1 x2 u_sym

    # Expansão de Taylor de -sin(x1 + θ_ref) + sin(θ_ref)
    termo_linear = -cos(θ_ref) * x1
    termo_quad   = (sin(θ_ref) / 2.0) * x1^2
    termo_cubico = (cos(θ_ref) / 6.0) * x1^3
    
    f_error = termo_linear + termo_quad + termo_cubico

    f = [ x2, f_error - β_ndim * x2 ]
    g = [0.0, 1.0]  
    h = [x1, x2]
# =========================================================

# =========================================================
# 4. SÍNTESE DO CONTROLADOR VIA DBC (no espaço adimensional)
    println("Sintetizando controlador DBC para o Pêndulo Simples...\n")

    K_ndim = otimizarDBC(
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
    # u_ndim = K_ndim · [x₁, x₂]
    # u_dim  = u_ndim · (m·L²·ω₀²)
    # u_dim  = K_ndim · (m·L²·ω₀²) · [x₁, x₂]
    #        = K_dim · [x₁, x₂]

    escala_u = m * L^2 * ω0^2    # = m·L²·ω₀² = m·g·L  (= I·ω₀²)

    K1 = K_ndim[1, 1] * escala_u
    K2 = K_ndim[1, 2] * escala_u

    @printf("\nGanhos no espaço adimensional : K_ndim = [%.4f  %.4f]\n", K_ndim[1,1], K_ndim[1,2])
    @printf("Ganhos físicos (torque [N·m]) : K₁   = %.4f   K₂ = %.4f\n", K1, K2)
    @printf("Lei de controle: u(x) = %.4f·θ + %.4f·ω\n\n", K1, K2)
# =========================================================

# =========================================================
# 6. STRUCT E DINÂMICAS PARA SIMULAÇÃO
    struct PenduloDBCParams
        L               ::Float64
        m               ::Float64
        β               ::Float64
        K1              ::Float64
        K2              ::Float64
        θ_ref           ::Float64
        estadosIniciais ::Vector{Float64}
        numEstados      ::Int
    end

    function PenduloDBCParams(; L, m, β, K1, K2, θ_ref, estadosIniciais)
        PenduloDBCParams(L, m, β, K1, K2, θ_ref, estadosIniciais, 2)
    end

    # Dinâmica exata em malha fechada 
    function pendulo_dbc_mf!(dx, x, p::PenduloDBCParams, t)
        I = p.m * p.L^2
        u = p.K1 * x[1] + p.K2 * x[2]

        dx[1] = x[2]
        dx[2] = (-p.m * 9.81 * p.L * sin(x[1]) - p.β * x[2] + u) / I
    end

    # Dinâmica exata em malha fechada 
    function pendulo_dbc_mf!(dx, x, p::PenduloDBCParams, t)
        I = p.m * p.L^2
        
        # 1. Mudança de coordenadas (Cálculo do Erro)
        x_tilde1 = x[1] - p.θ_ref
        x_tilde2 = x[2] - 0.0 # A velocidade alvo é sempre zero
        
        # 2. Esforço dinâmico realimentado (SOS)
        u_dinamico = p.K1 * x_tilde1 + p.K2 * x_tilde2
        
        # 3. Esforço de regime (vencer a gravidade no alvo)
        # Torque físico constante necessário: m * g * L * sin(θ_ref)
        u_ss = p.m * 9.81 * p.L * sin(p.θ_ref)
        
        # 4. Controle total aplicado no motor
        u = u_dinamico + u_ss

        dx[1] = x[2]
        # Aplica a física EXATA (sin(x[1])) usando o controle total
        dx[2] = (-p.m * 9.81 * p.L * sin(x[1]) - p.β * x[2] + u) / I
    end

    # Dinâmica em malha aberta
    function pendulo_aberta!(dx, x, p::PenduloDBCParams, t)
        I = p.m * p.L^2
        dx[1] = x[2]
        dx[2] = (-p.m * 9.81 * p.L * sin(x[1]) - p.β * x[2]) / I
    end
# =========================================================

# =========================================================
# 7. SIMULAÇÃO
    params = PenduloDBCParams(L  = L,
                              m  = m,
                              β  = β,
                              K1 = K1,
                              K2 = K2,
                              θ_ref = θ_ref,
                              estadosIniciais = x0)

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
                                  referencia = θ_ref,
                                  idx_estado = 1,
                                  banda      = 0.05)

    m_omega = analisarPerformance(sol_fechada;
                                  referencia = θ_ref,
                                  idx_estado = 2,
                                  banda      = 0.05)

    imprimirRelatorio(m_theta, nome = "Pêndulo DBC — Ângulo θ")
    imprimirRelatorio(m_omega, nome = "Pêndulo DBC — Velocidade angular ω")
# =========================================================

# =========================================================
# 10. COMPARAÇÃO COM PID
    include("../../Controladores/PID.jl")
    include("../../Abstrações/MalhaFechada.jl")

    pendulo_planta = Pendulo(
        comprimento     = L,
        massa           = m,
        coefAtrito      = β,
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