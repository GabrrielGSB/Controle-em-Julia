# Controle/Exemplos/MEMS/MEMS_DBC.jl
#
# Estabilização do Sistema Microeletromecânico (MEMS) via
# Dissipativity-Based Conditions (DBC) com otimização SOS.
#
# Referência: "Necessary and Sufficient Dissipativity-Based Conditions
#              for Feedback Stabilization" — Diego de S. Madeira
#
# Fluxo:
#   1. Obtém o ganho K via otimizarDBC()
#   2. Fecha a malha com u = K * x3
#   3. Simula e compara malha aberta × malha fechada
#   4. Analisa performance do controlador encontrado

# =========================================================
# INCLUDES
    include("../SistemasLTI/MEMS.jl")
    include("../Ferramentas/ResoluçãoEDO.jl")
    include("../Ferramentas/Visualização.jl")
    include("../Ferramentas/AnálisePerformance.jl")
    include("../Ferramentas/DBC.jl")
# =========================================================


# =========================================================
# 1. PARÂMETROS DA PLANTA
# =========================================================
r_param = 1.0
ξ       = 0.1
q̄       = 1.0
x̄       = 0.5

# Condições iniciais: pequena perturbação em torno da origem
x0     = [0.1, -0.05, 0.2]
tspan  = (0.0, 15.0)


# =========================================================
# 2. SÍNTESE DO CONTROLADOR VIA DBC (SOS)
# =========================================================
# O controlador tem a forma u = K * h(x), onde h(x) = x3
# (a carga elétrica é a saída disponível para realimentação)

@polyvar x1 x2 x3 u_sym

f_mems = [
    x2,
    -x1 - 2*ξ*x2 + (2*q̄/3)*x3 + (1/3)*x3^2,
    (q̄/r_param)*x1 + (1/r_param)*x1*x3 + ((x̄ - 1.0)/r_param)*x3
]
g_mems = [0.0, 0.0, 2.0/(3.0*r_param)]   # Coluna de entrada — atua apenas em dx3
h_mems = [x3]                              # Saída fictícia: y = x3

K_sos = otimizarDBC(
    f_mems, g_mems, h_mems,
    [x1, x2, x3], u_sym;
    grau_V    = 1:2,
    grau_T    = 1:4,
    epsilon_R = 1e-4,
    beta      = 1e-8,
    verbose   = true
)

# K_sos é uma Matrix{Float64} (1×1) — extrai o escalar
K_dbc = K_sos[1, 1]
println("Ganho DBC encontrado: K = $K_dbc\n")


# =========================================================
# 3. SIMULAÇÃO — MALHA ABERTA vs. MALHA FECHADA
# =========================================================

# --- Malha Aberta (K = 0) ---
params_aberta = MEMSParams(
    q_bar   = q̄,
    xi      = ξ,
    r_param = r_param,
    x_bar   = x̄,
    K       = 0.0
)

sol_aberta = resolverSistema(MEMS!, x0, tspan, params_aberta; resolucao=0.05)

# --- Malha Fechada (K calculado via DBC) ---
params_fechada = MEMSParams(
    q_bar   = q̄,
    xi      = ξ,
    r_param = r_param,
    x_bar   = x̄,
    K       = K_dbc
)

sol_fechada = resolverSistema(MEMS!, x0, tspan, params_fechada; resolucao=0.05)


# =========================================================
# 4. VISUALIZAÇÃO
# =========================================================

# 4.1 Evolução temporal: estado de deflexão (x1) — malha aberta vs fechada
plotarNoTempo(
    [sol_aberta, sol_fechada];
    nomes   = ["Malha Aberta", "DBC (SOS)"],
    titulo  = "MEMS — Deflexão x₁: Malha Aberta vs. DBC",
    estados = 1:3
)

# 4.2 Retrato de fase: deflexão × velocidade (x1 × x2)
plotarRetratoFase(
    [sol_aberta, sol_fechada];
    estados = (1, 2),
    nomes   = ["Malha Aberta", "DBC (SOS)"],
    titulo  = "MEMS — Retrato de Fase (Deflexão × Velocidade)"
)

# 4.3 Retrato de fase: deflexão × carga (x1 × x3)
plotarRetratoFase(
    [sol_aberta, sol_fechada];
    estados = (1, 3),
    nomes   = ["Malha Aberta", "DBC (SOS)"],
    titulo  = "MEMS — Retrato de Fase (Deflexão × Carga)"
)

# 4.4 Campo vetorial da malha fechada no plano (x1 × x2)
#     (fixa x3 = 0 para visualização 2D do plano principal)
function mems_plano_x1x2!(dx, x, p, t)
    x_3d = [x[1], x[2], 0.0]
    dx_3d = zeros(3)
    MEMS!(dx_3d, x_3d, p, t)
    dx[1] = dx_3d[1]
    dx[2] = dx_3d[2]
end

plotarRetratoFaseCompleto(
    mems_plano_x1x2!, params_fechada;
    limiteEixo         = 0.5,
    densidadeSetas     = 16,
    raiosIniciais      = [0.1, 0.25, 0.4],
    trajetoriasPorAnel = 6,
    tempoMaximo        = 15.0,
    titulo             = "MEMS — Campo Vetorial (Malha Fechada, x₃=0)"
)


# =========================================================
# 5. ANÁLISE DE PERFORMANCE
# =========================================================

# O MEMS é regulador: referência = 0 (todos os estados → 0)

m_x1 = analisarPerformance(sol_fechada;
    referencia  = 0.0,
    idx_estado  = 1,
    banda       = 0.05   # ±5% pois x0 é pequeno
)

m_x2 = analisarPerformance(sol_fechada;
    referencia = 0.0,
    idx_estado = 2,
    banda      = 0.05
)

m_x3 = analisarPerformance(sol_fechada;
    referencia = 0.0,
    idx_estado = 3,
    banda      = 0.05
)

imprimirRelatorio(m_x1, nome = "MEMS DBC — Deflexão x₁")
imprimirRelatorio(m_x2, nome = "MEMS DBC — Velocidade x₂")
imprimirRelatorio(m_x3, nome = "MEMS DBC — Carga x₃")

# Comparação direta com o controlador analítico do artigo
S_art = -(x̄ - 1.0) / 3.0
R_art = 0.1
K_art = -(1.0 / R_art) * S_art

println("\n─── Comparação de Ganhos ─────────────────────────────")
@printf("  K  analítico (artigo) = %+.4f\n", K_art)
@printf("  K  DBC (SOS)          = %+.4f\n", K_dbc)
println("──────────────────────────────────────────────────────\n")

params_art = MEMSParams(
    q_bar   = q̄,
    xi      = ξ,
    r_param = r_param,
    x_bar   = x̄,
    K       = K_art
)
sol_art = resolverSistema(MEMS!, x0, tspan, params_art; resolucao=0.05)

compararControladores(
    [sol_art, sol_fechada],
    ["Analítico (artigo)", "DBC (SOS)"];
    referencia = 0.0,
    idx_estado = 1,       # analisa deflexão x1
    banda      = 0.05
)