using ForwardDiff
using LinearAlgebra
using Printf

# Importação das ferramentas e abstrações
include("../../Controladores/RE.jl")
include("../../Abstrações/MalhaFechada.jl")
include("../../SistemasLTI/PênduloInvertido.jl")
include("../../Ferramentas/ResoluçãoEDO.jl")
include("../../Ferramentas/Visualização.jl")
include("../../Ferramentas/projetarPolos.jl")
include("../../Ferramentas/LinearizarSistema.jl")
include("../../Ferramentas/AnálisePerformance.jl")
include("../../Algoritmos/Bass_Gura.jl")

# =========================================================================
# CONFIGURAÇÃO DA PLANTA E CONTROLE
    # Instancia a planta (parâmetros padrão)
    planta = PenduloInvertido(M=1.0, m=0.5, l=0.5, g=9.81, b_c=40, b_p=0.1)

    # Ponto de Operação (Equilíbrio instável: haste para cima)
    x_eq = [0.0, 0.0, 0.0, 0.0]
    u_eq = [0.0]

    # Linearização automática
    A, B_mat, _, _ = linearizar(planta, x_eq, u_eq)
    B = B_mat[:, 1]

    # Projeto de Ganhos via Bass-Gura
    polos_desejados = projetarPolos(4, 10.0, 3.0) 
    K = bass_gura(A, B, polos_desejados)
    K = K'
# =========================================================================

# =========================================================================
# MONTAGEM DAS MALHAS ABERTA E FECHADA
    ref = [1.0, 0.0, 0.0, 0.0]  

    # A) Malha Aberta (Ganho K = 0) -> Ação de controle será sempre nula
    ctrl_aberto = RealimentacaoEstado(zeros(1, 4))
    sys_aberto  = conectar(planta, ctrl_aberto, ref)

    # B) Malha Fechada (Ganho K calculado por Bass-Gura)
    ctrl_fechado = RealimentacaoEstado(K)
    sys_fechado  = conectar(planta, ctrl_fechado, ref)

# =========================================================================

# =========================================================================
# SIMULAÇÕES
    # Condição inicial: carrinho no zero, haste levemente inclinada (10 graus = 0.17 rad)
    x0_fisico = [0.0, 0.0, 0.17, 0.0]
    x0_total  = condicoesIniciais(sys_fechado, x0_fisico) # A mesma dimensão para ambos

    tspan = (0.0, 10.0)

    println("\nSimulando Malha Aberta...")
    sol_aberta = resolverSistema(sys_aberto, x0_total, tspan)

    println("Simulando Malha Fechada...")
    sol_fechada = resolverSistema(sys_fechado, x0_total, tspan)
# =========================================================================

# =========================================================================
# VISUALIZAÇÃO E ANÁLISE
    println("\n--- Relatório de Estabilização ---")
    println("Polos de malha aberta: ", round.(eigvals(A), digits=2))
    println("Polos de malha fechada reais: ", round.(eigvals(A - B*K), digits=2))

    # 4.1 Gráficos da Malha Aberta
    plotarNoTempo(sol_aberta, titulo="Malha Aberta - Pêndulo Caindo", estados=(1, 3))
    plotarRetratoFase(sol_aberta, estados=(1,3))
    # plotarRetratoFaseCompleto(sol_aberta, titulo="Retrato de Fase - Malha Aberta")

    # 4.2 Gráficos da Malha Fechada
    plotarNoTempo(sol_fechada, titulo="Malha Fechada - Estabilização (Bass-Gura)", estados=(1, 3))
    plotarRetratoFase(sol_fechada, estados=(1,3))
    # plotarRetratoFaseCompleto(sol_fechada, titulo="Retrato de Fase - Malha Fechada")

    metricas = analisarPerformance(sol_fechada, referencia=0, idx_estado=3)
    imprimirRelatorio(metricas)
# =========================================================================