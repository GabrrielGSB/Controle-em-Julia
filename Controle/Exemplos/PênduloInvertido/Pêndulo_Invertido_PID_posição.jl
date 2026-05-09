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
    planta = PenduloInvertido(M = 1.0, m = 0.1, l = 0.5, g = 9.81, b = 0.1)

    pid_posicao = PID(Kp = 5.0, 
                      Ki = 0.1, 
                      Kd = 2.0)
# =========================================================

# =========================================================
# SIMULAÇÃO
    sys = conectar(planta      = planta, 
                   controlador = pid_posicao, 
                   referencia  = 1.5, 
                   idx_saida   = [1, 2])  
    
    x0  = condicoesIniciais(sys, [0.0, 0.0, 0.0, 0.0]) 
    
    t_simu = 15.0  
    sol = resolverSistema(sys, x0, (0.0, t_simu),
                          resolucao=0.01,
                          salvar_controle=false)
# =========================================================

# =========================================================
# VISUALIZAÇÃO E ANÁLISE
    plotarNoTempo(sol, titulo="PID no Pêndulo Invertido (Controle de POSIÇÃO)", 
                  estados=(1, 3), # Focamos na Posição (1)
                  mostrar_controle=false)

    # Analisamos a performance da Posição (estado 1)
    m = analisarPerformance(sol, referencia=1.5, idx_estado=1)
    imprimirRelatorio(m)
# =========================================================