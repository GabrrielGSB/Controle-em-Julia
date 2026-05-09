#=
SIMULAÇÃO EM MALHA ABERTA: PÊNDULO INVERTIDO NO CARRINHO
-----------------------------------------------------------------
Objetivo: Observar a instabilidade natural do sistema quando 
nenhuma força de controle é aplicada. O pêndulo deve cair e 
movimentar o carrinho por reação.
=#

# 1. Importação dos Módulos do Framework
include("../../SistemasLTI/PênduloInvertido.jl")
include("../../Ferramentas/ResoluçãoEDO.jl")
include("../../Ferramentas/Visualização.jl")

println("Iniciando simulação de Malha Aberta do Cart-Pole...")

#=========================================================================
2. DEFINIÇÃO DA PLANTA E CONDIÇÕES INICIAIS
=========================================================================#

# Vamos iniciar o carrinho parado na origem, mas com o pêndulo levemente 
# inclinado (0.1 radianos) para induzir a queda.
x0_instavel = [
    0.0,  # x[1] : Posição do carrinho (m)
    0.0,  # x[2] : Velocidade do carrinho (m/s)
    0.1,  # x[3] : Ângulo do pêndulo (rad) -> Quase em pé
    0.0   # x[4] : Velocidade angular (rad/s)
]

# Instanciando a struct da planta com os valores padrão
planta_cartpole = PenduloInvertido(estadosIniciais = x0_instavel)

# Tempo de simulação curto, pois a queda é rápida
tspan = (0.0, 5.0)

#=========================================================================
3. RESOLUÇÃO DA EDO
=========================================================================#

# O resolverSistema vai integrar a dinâmica pendulo_invertido! que 
# definimos no arquivo da planta.
solucao_aberta = resolverSistema(
    pendulo_invertido!, 
    planta_cartpole.estadosIniciais, 
    tspan, 
    planta_cartpole; 
    resolucao = 0.02
)

#=========================================================================
4. VISUALIZAÇÃO DOS RESULTADOS (PlotlyJS)
=========================================================================#

# Vamos observar a Posição do Carrinho (estado 1) e o Ângulo (estado 3)
plotarNoTempo(
    solucao_aberta;
    titulo = "Cart-Pole: Queda Livre (Malha Aberta)",
    estados = (1, 3)
)

plotarRetratoFase(
    solucao_aberta;
    estados = (3, 4),
    titulo = "Cart-Pole — Retrato de Fase (Ângulo vs. Velocidade Angular)"
)

plotarRetratoFase(
    solucao_aberta;
    estados = (1, 3), # x_1 (Posição) e x_3 (Ângulo)
    titulo = "Cart-Pole — Dinâmica Real (Posição Carro vs. Ângulo)"
)


#=========================================================================
5. CAMPO VETORIAL COMPLETO (Retrato de Fase 2D)
=========================================================================#
println("Gerando o Campo Vetorial Completo (isso pode levar alguns segundos)...")

# Função Wrapper: Converte o plano 2D (θ, ω) para o sistema 4D completo
function cartpole_plano_x3x4!(dx, x, p, t)
    # x[1] recebido do plot é o ângulo (antigo x[3])
    # x[2] recebido do plot é a velocidade angular (antigo x[4])
    
    # Monta o vetor 4D "congelando" a posição e velocidade do carrinho em zero
    x_4d = [0.0, 0.0, x[1], x[2]]
    dx_4d = zeros(4)
    
    # Chama a dinâmica real que definimos na struct
    pendulo_invertido!(dx_4d, x_4d, p, t; u=0.0)
    
    # Extrai apenas as derivadas rotacionais para o gráfico 2D
    dx[1] = dx_4d[3]
    dx[2] = dx_4d[4]
end

# Gerando o plot do campo vetorial
plotarRetratoFaseCompleto(
    cartpole_plano_x3x4!, 
    planta_cartpole; 
    limiteEixo = 3.14,            # Vai até ~3.14 radianos (pêndulo caído para baixo)
    densidadeSetas = 16,          # Reduzido um pouco para não poluir a tela
    raiosIniciais = [0.1, 1.0, 2.5], # Lança trajetórias perto do topo, meio e fundo
    trajetoriasPorAnel = 4, 
    tempoMaximo = 4.0, 
    titulo = "Cart-Pole — Campo Vetorial Completo (Ângulo vs. Vel. Angular)"
)

println("Simulação e visualizações concluídas com sucesso!")


