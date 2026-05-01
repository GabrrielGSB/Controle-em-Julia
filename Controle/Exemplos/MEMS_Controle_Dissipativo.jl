# Implementado de acordo com artigo: Necessary and Sufficient Dissipativity-Based Conditions for Feedback Stabilization
# por Diego de S. Madeira
using Plots

include("../SistemasLTI/MEMS.jl")
include("../Ferramentas/ResoluçãoEDO.jl")

# ==========================================
# 1. Definição dos Parâmetros da Planta
# ==========================================
r = 1.0       
ξ = 0.1      
q̄ = 1.0       
x̄ = 0.5       

# ==========================================
# 2. Projeto do Controlador (Dissipatividade)
# ==========================================
S = -(x̄ - 1.0) / 3.0  # Conforme Equação (78) do artigo:
R = 0.1               # Escolhemos um R > 0 pequeno:                 
K = -(1.0 / R) * S    # O ganho do controlador estabilizante K = -R^(-1) * S^T  

parametros = MEMSParams(q_bar   = q̄, 
                        xi      = ξ, 
                        r_param = r, 
                        x_bar   = x̄,
                        K       = K)

# ==========================================
# 4. Simulação
# ==========================================
# Condições Iniciais: perturbamos o sistema tirando-o da origem
x0 = [0.1, -0.05, 0.2]  # [Deflexão, Velocidade, Carga]
tspan = (0.0, 15.0)


solucao = resolverSistema(MEMS!, x0, tspan, parametros)

# # Criando o problema de Equações Diferenciais Ordinárias (ODE)
# prob = ODEProblem(MEMS!, x0, tspan)

# # Resolvendo numericamente com o algoritmo Tsit5 (Runge-Kutta)
# sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# ==========================================
# 5. Plotagem dos Resultados
# ==========================================
plot(solucao, 
     labels=["x1 (Deflexão)" "x2 (Velocidade)" "x3 (Carga)"],
     title="Estabilização do MEMS via QSR-Dissipatividade",
     xlabel="Tempo (s)", 
     ylabel="Amplitude dos Estados",
     linewidth=2,
     grid=true)