include("FuncSimuEDO.jl")
include("SistemasLTI.jl")

paramSistema = PenduloParams(gravidade   = 9.81, 
                             comprimento = 1.0, 
                             massa       = 1.0, 
                             coefAtrito  = 0.5,
                             estadosIniciais = (θ = 2.0, ω = 0.0)) 

# paramSistema = MassaMolaParams(coefElastico    = 40.0, 
#                                coefAtrito      = 0.8, 
#                                massa           = 1.0, 
#                                estadosIniciais = (x = 0.5, v = 0.0))                            

simularSistemaMalhaAberta(sistema         = pendulo!, 
                          estadosIniciais = paramSistema.estadosIniciais, 
                          parametros      = paramSistema, 
                          tempoSimulacao  = (0.0, 10.0),
                          visualizacao    = true,
                          animar          = true,
                          animarFPS       = 24)
