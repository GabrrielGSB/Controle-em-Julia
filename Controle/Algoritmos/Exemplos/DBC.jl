using JuMP
using SumOfSquares
using DynamicPolynomials
using LinearAlgebra
using COSMO
using Clarabel
# import SCS # Solver SDP (equivalente open-source ao SeDuMi)

# ======================================================================== #
#                   IMPLEMENTAÇÃO SOS - RECURRENT DBIs                     #
#                FUNCIONA MUITO BEM: DELTA_k+1 > DELTA_k                   #
#                     ESTABILIZAÇÃO GLOBAL, SEM EID                        #
#                EXEMPLO 1: MADEIRA E MACHADO MICNON 2024                  #
# ======================================================================== #

println("Iniciando a Otimização SOS...")
# --------------- DEFINIÇÃO DA PLANTA / LEI CONTROLE ------------------- #
# 1. Declaração das Variáveis Polinomiais
@polyvar x1 x2 u
x = [x1, x2]

# Dinâmica da planta (Exemplo 2 do artigo Valmorbida and Papachristodoulou)
f = [2*x1^3 + (x1^2)*x2 - 6*x1*x2^2 + 5*(x2^3), 0]
g = [0, 1]

# Saída fictícia h(x)
h_poly = x1^3 + x1^2*x2 + x1*x2^2 + x2^3
y = [h_poly] # Mantido como vetor de 1 elemento para operações de matriz

# ------------------------- PRIMEIRA TENTATIVA ------------------------- #
# Inicialização do modelo SOS (STEP 1)
# model0 = SOSModel(SCS.Optimizer)
# model0 = SOSModel(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "eps_abs" => 1e-5))
model0 = SOSModel(Clarabel.Optimizer)
set_silent(model0)

# Função candidata de Lyapunov V(x) e função T(x)
Z_v = monomials(x, 1:2)
@variable(model0, V0, SOSPoly(Z_v))
@variable(model0, T0, SOSPoly(Z_v))

# Matrizes reais constantes Q, S, R (neste exemplo, 1x1)
@variable(model0, Q0[1:1, 1:1], Symmetric)
@variable(model0, S0[1:1, 1:1])
@variable(model0, R0[1:1, 1:1], Symmetric)

epsilon_R = 1e-4
@constraint(model0, R0[1,1] >= epsilon_R)

# Restrições de positividade estrita (margens)
beta = 1e-8
@constraint(model0, V0 - beta*(x1^2 + x2^2)   in SOSCone())
@constraint(model0, T0 - beta*(x1^2 + x2^2)^2 in SOSCone())

# Gradiente de V0
gradV0 = differentiate(V0, x)

# Desigualdade de dissipatividade estrita: -Vdot - T + r(y,u) >= 0
Vdot0 = dot(gradV0, f + g*u)
r_yu0 = dot(y, Q0 * y) + 2 * dot(y, S0 * [u]) + dot([u], R0 * [u])
diss_expr0 = -Vdot0 - T0 + r_yu0

@constraint(model0, diss_expr0 in SOSCone())

# Resolvendo o SOSP inicial
optimize!(model0)

if termination_status(model0) != MOI.OPTIMAL && termination_status(model0) != MOI.ALMOST_OPTIMAL
    error("O solver não encontrou uma solução inicial factível na Iteração 0.")
end

ITERACAO = 0
println("ITERACAO = ", ITERACAO)

# Extraindo os valores numéricos da iteração 0
Q0_val = value.(Q0)
S0_val = value.(S0)
R0_val = value.(R0)

# Cálculo do Delta Inicial
Delta_0 = S0_val*inv(R0_val)*S0_val' - Q0_val
EIG_DELTA_0 = eigvals(Delta_0)[1]
println("EIG_DELTA_0 = ", EIG_DELTA_0)

# -------------------------- PARTE ITERATIVA --------------------------- #
imax = 100
i = 1

# Variáveis para guardar o resultado final
Q_final, S_final, R_final = Q0_val, S0_val, R0_val

while i < imax
    global Q0_val, S0_val, R0_val, Q_final, S_final, R_final, i, imax
    
    # model = SOSModel(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "eps_abs" => 1e-5))
    model = SOSModel(Clarabel.Optimizer)
    set_silent(model)
    
    # Variáveis da iteração k+1
    @variable(model, V, SOSPoly(Z_v))
    @variable(model, T, SOSPoly(Z_v))
    
    @variable(model, Q[1:1, 1:1], Symmetric)
    @variable(model, S[1:1, 1:1])
    @variable(model, R[1:1, 1:1], Symmetric)

    @constraint(model, R[1,1] >= epsilon_R)
    
    # Restrições de positividade
    @constraint(model, V - beta*(x1^2 + x2^2) in SOSCone())
    @constraint(model, T - beta*(x1^2 + x2^2)^2 in SOSCone())
    
    # Gradiente e dissipatividade
    gradV = differentiate(V, x)
    Vdot = dot(gradV, f + g*u)
    r_yu = dot(y, Q * y) + 2 * dot(y, S * [u]) + dot([u], R * [u])
    diss_expr = -Vdot - T + r_yu
    
    @constraint(model, diss_expr in SOSCone())
    
    # Restrições Recorrentes (LMIs) para garantir o crescimento de Delta_c
    M = S*inv(R0_val)*S0_val' + S0_val*inv(R0_val)*S' - 2*S0_val*inv(R0_val)*S0_val' + Q0_val - Q
    
    # Em JuMP, restrições LMI usam PSDCone (Positive Semidefinite Cone)
    @constraint(model, Symmetric(M) in PSDCone())
    @constraint(model, Symmetric(R0_val - R) in PSDCone())
    
    # Resolver
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL && termination_status(model) != MOI.ALMOST_OPTIMAL
        println("Aviso: Solver falhou na iteração $i. Interrompendo loop.")
        break
    end
    
    # Extrair novos valores
    Q_val = value.(Q)
    S_val = value.(S)
    R_val = value.(R)
    
    Delta_i = S_val * inv(R_val) * S_val' - Q_val
    EIG_DELTA_i = eigvals(Delta_i)[1]
    
    println("ITERACAO = ", i)
    println("EIG_DELTA_i = ", round(EIG_DELTA_i, digits=5))
    
    # Condição de Parada: O problema foi estabilizado!
    if EIG_DELTA_i >= 0
        Q_final = Q_val
        S_final = S_val
        R_final = R_val
        println("Solução estabilizante encontrada!")
        break
    end
    
    # Atualiza as matrizes para a próxima iteração
    Q0_val = Q_val
    S0_val = S_val
    R0_val = R_val
    
    i += 1
end

# -------------------------- RESULTADOS E CONTROLE --------------------- #
K = -inv(R_final) * S_final'
println("Ganho do Controlador (K): ", round(K[1,1], digits=4))

# A lei de controle em malha fechada seria uk = K[1,1] * h_poly


using OrdinaryDiffEq
using PlotlyJS
using LinearAlgebra

# 1. Definindo a Dinâmica do Sistema
function sistema_polinomial!(dx, x, p, t)
    x1, x2 = x[1], x[2]
    K, controle_ativo = p # p agora será uma Tupla, evitando o erro de tipagem
    
    # Dinâmica em malha aberta
    f1 = 2*x1^3 + (x1^2)*x2 - 6*x1*x2^2 + 5*x2^3
    f2 = 0.0
    
    # Saída fictícia h(x)
    h = x1^3 + x1^2*x2 + x1*x2^2 + x2^3
    
    # Aplicação do controle u(x)
    u = controle_ativo ? K * h : 0.0
    
    # Derivadas
    dx[1] = f1
    dx[2] = f2 + u
end

# 2. Parâmetros e Condições Iniciais
K_calculado = -4.7834 
# Usando TUPLAS (parênteses) para misturar Float e Bool sem problemas de conversão
p_aberta = (K_calculado, false) 
p_fechada = (K_calculado, true)

x0 = [0.5, -0.5] 
tspan = (0.0, 5.0)

# 3. Resolvendo as Equações Diferenciais
# Malha fechada (Estável)
prob_fechada = ODEProblem(sistema_polinomial!, x0, tspan, p_fechada)
sol_fechada = solve(prob_fechada, Tsit5(), reltol=1e-6, abstol=1e-6)

# Malha aberta (Instável - Callback para abortar antes do infinito)
condicao_explosao(u, t, integrador) = norm(u) > 100.0
terminar_simulacao!(integrador) = terminate!(integrador)
cb = DiscreteCallback(condicao_explosao, terminar_simulacao!)

prob_aberta = ODEProblem(sistema_polinomial!, x0, tspan, p_aberta)
sol_aberta = solve(prob_aberta, Tsit5(), callback=cb, reltol=1e-6, abstol=1e-6)


# ======================================================================== #
#                        VISUALIZAÇÃO COM PLOTLYJS                         #
# ======================================================================== #

# Extraindo os dados das soluções
t_f = sol_fechada.t
x1_f = [u[1] for u in sol_fechada.u]
x2_f = [u[2] for u in sol_fechada.u]

t_a = sol_aberta.t
x1_a = [u[1] for u in sol_aberta.u]
x2_a = [u[2] for u in sol_aberta.u]

# Criando as linhas (Traces) para o Gráfico 1: Tempo
trace_x1_f = scatter(x=t_f, y=x1_f, mode="lines", name="x1 (Com Controle)", line_color="red")
trace_x2_f = scatter(x=t_f, y=x2_f, mode="lines", name="x2 (Com Controle)", line_color="blue")

# Criando as linhas (Traces) para o Gráfico 2: Plano de Fase
trace_fase_aberta = scatter(x=x1_a, y=x2_a, mode="lines", name="Sem Controle (Diverge)", 
                            line=attr(dash="dash", color="orange", width=2))
trace_fase_fechada = scatter(x=x1_f, y=x2_f, mode="lines", name="Com Controle (Estável)", 
                             line=attr(color="green", width=2))

# Marcadores especiais para o Início e para a Origem
trace_inicio = scatter(x=[x0[1]], y=[x0[2]], mode="markers", name="Início (x0)", 
                       marker=attr(color="black", size=10))
trace_origem = scatter(x=[0], y=[0], mode="markers", name="Origem (0,0)", 
                       marker=attr(color="purple", size=12, symbol="star"))

# Configurando o Layout com Subplots (1 linha, 2 colunas)
fig = make_subplots(
    rows=1, cols=2, 
    subplot_titles=["<b>Estados x(t) - Malha Fechada</b>" "<b>Retrato de Fase (x1 vs x2)</b>"]
)

# Adicionando os Traces aos painéis respectivos
add_trace!(fig, trace_x1_f, row=1, col=1)
add_trace!(fig, trace_x2_f, row=1, col=1)

add_trace!(fig, trace_fase_aberta, row=1, col=2)
add_trace!(fig, trace_fase_fechada, row=1, col=2)
add_trace!(fig, trace_inicio, row=1, col=2)
add_trace!(fig, trace_origem, row=1, col=2)

# Melhorando a aparência geral e nomeando os eixos
relayout!(fig, 
    title_text="Dashboard de Controle: Estabilização de Sistema Polinomial",
    title_font_size=20,
    xaxis_title="Tempo (s)",
    yaxis_title="Amplitude",
    xaxis2_title="x1",
    yaxis2_title="x2",
    height=550, 
    width=1100,
    template="plotly_white",
    hovermode="closest"
)

# Exibe o gráfico interativo
display(fig)