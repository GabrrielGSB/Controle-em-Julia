include("FuncSimuEDO.jl")
include("SistemasLTI.jl")

simularSistemaMalhaAberta(sistema         = pendulo!, 
                          estadosIniciais = [1, 0], 
                          parametros      = (g = 9.81, L = 1, m = 1, b = 1), 
                          tempoSimulacao  = (0.0, 30.0))