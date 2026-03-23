include("FuncSimuEDO.jl")
include("SistemasLTI.jl")

paramSistema = PenduloParams(gravidade   = 9.81, 
                             comprimento = 1.0, 
                             massa       = 1.0, 
                             coefAtrito  = 0.5,
                             estadosIniciais = (θ = 1.0, ω = 0.0)) 

simularSistemaMalhaAberta(sistema         = pendulo!, 
                          estadosIniciais = paramSistema.estadosIniciais, 
                          parametros      = paramSistema, 
                          tempoSimulacao  = (0.0, 10.0),
                          visualizacao    = true,
                          animar          = true)
