# =========================================================
# INCLUDES
    include("../../SistemasLTI/PênduloInvertido.jl")
    include("../../Controladores/PID.jl")
    include("../../Abstrações/MalhaFechada.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

# =========================================================
# PARÂMETROS FÍSICOS E DE CONTROLE
    planta = PenduloInvertido(M = 1.0, 
                              m = 0.1, 
                              l = 0.5, 
                              g = 9.81, 
                              b_c = 0.1,
                              b_p = 0.1)

    # Ganhos negativos pois a força do carrinho deve contrariar a queda da haste
    pid = PID(Kp = -50.0,
              Ki = -1.0,
              Kd = -10.0)
# =========================================================

# =========================================================
# MONTAGEM DA MALHA FECHADA
    sys = conectar(planta, pid, 0.0, [3, 4])  
    
    # Estados iniciais: [x, v, theta, omega] -> Iniciando com 0.1 rad de inclinação
    x0  = condicoesIniciais(sys, [0.0, 0.0, 0.1, 0.0]) 
# =========================================================

# =========================================================
# SIMULAÇÃO 
    t_simu = 15.0  
    sol = resolverSistema(sys, x0, (0.0, t_simu),
                          resolucao=0.01,
                          salvar_controle=false)
# =========================================================

# =========================================================
# VISUALIZAÇÃO E ANÁLISE
    plotarNoTempo(sol, titulo="PID no Pêndulo Invertido (Controle de Ângulo)", 
                  estados=(1, 3), # Plota Posição (1) e Ângulo (3) simultaneamente
                  mostrar_controle=false)

    # m = analisarPerformance(sol, referencia=0.0, idx_estado=3); imprimirRelatorio(m)

    # gerarAnimacao(sol, planta, 30, mostrar_controle=false)
# =========================================================