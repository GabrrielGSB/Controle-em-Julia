# =========================================================
# RE-APROVEITAMENTO DE CÓDIGO
    include("../SistemasLTI/Pêndulo.jl")
    include("../Ferramentas/ResoluçãoEDO.jl")
    include("../Ferramentas/Visualização.jl")
# =========================================================

# =========================================================
# PARÂMETROS FÍSICOS
    parametros = PenduloParams(
        comprimento     = 1.0,   # L [m]
        massa           = 1.0,   # m [kg]
        coefAtrito      = 0.6,   # b [N.m.s/rad]
        estadosIniciais = [π/4, 0.0]  # [θ₀ (rad), ω₀ (rad/s)]
    )
# =========================================================

# =========================================================
# SIMULAÇÃO (malha aberta)
    intervaloTempo = (0.0, 10.0)

    solucao = resolverSistema(pendulo!,
                              parametros.estadosIniciais,
                              intervaloTempo,
                              parametros)
# =========================================================

# =========================================================
# VISUALIZAÇÃO 
    plotarNoTempo(solucao,
                  titulo  = "Pêndulo Simples — Evolução Temporal",
                  estados = 1:2)

    plotarRetratoFase(solucao,
                      estados = (1, 2),
                      titulo  = "Pêndulo Simples — Retrato de Fase (θ × ω)")

    plotarRetratoFaseCompleto(pendulo!, parametros,
                              limiteEixo         = 3.0,
                              densidadeSetas     = 18,
                              raiosIniciais      = [0.5, 1.5, 2.5],
                              trajetoriasPorAnel = 1,
                              tempoMaximo        = 20.0,
                              titulo             = "Retrato de Fase - Campo Vetorial")
# =========================================================

# =========================================================
# ANIMAÇÃO

   gerarAnimacao(solucao, parametros, 30)    

# =========================================================