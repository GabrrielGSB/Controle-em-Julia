# =========================================================
# INCLUDES
    include("../../SistemasLTI/Pêndulo.jl")
    include("../../controladores/PID.jl")
    include("../../Abstrações/MalhaFechada.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

# =========================================================
# PARÂMETROS FÍSICOS E DE CONTROLE
    pendulo = Pendulo(comprimento     = 1.0,
                      massa           = 10.0,
                      coefAtrito      = 0.1,
                      estadosIniciais = [1.0, 0.0])

    pid = PID(Kp = 200.0,
              Ki = 20.0,
              Kd = 10.0)
# =========================================================

# =========================================================
# SIMULAÇÃO
    sys = conectar(pendulo, pid, π)      
    x0  = condicoesIniciais(sys, [1.0, 0.0])   
    sol = resolverSistema(sys, x0, (0.0, 10.0),
                          resolucao=0.01,
                          salvar_controle=false)
# =========================================================

# =========================================================
# VISUALIZAÇÃO E ANÁLISE
    plotarNoTempo(sol, titulo="PID no Pêndulo Simples", 
                  estados=(1), 
                  mostrar_controle=false)

    # m = analisarPerformance(sol, referencia=π, idx_estado=1); imprimirRelatorio(m)

    # gerarAnimacao(sol, pendulo, 30, mostrar_controle=false)
# =========================================================