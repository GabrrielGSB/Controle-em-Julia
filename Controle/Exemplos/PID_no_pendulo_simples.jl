# include("../SistemasLTI/Pêndulo.jl")
include("../Funções Auxiliares/Visualização.jl")
include("../Funções Auxiliares/ResoluçãoEDO.jl")
include("../SistemasLTI/Pêndulo.jl")



parametros_fisicos = PenduloParams(gravidade    = 9.81, 
                                   comprimento  = 1.0, 
                                   massa        = 1.0, 
                                   coefAtrito   = 0.1, 
                                   estadosIniciais = [1, 0.0, 0.0])


parametros_controlador_PID = PIDparams{PenduloParams}(parametros_fisicos, 
                                                      200.0,  # Kp
                                                      20.0,   # Ki
                                                      10.0,   # Kd
                                                      π)

intervaloTempo = (0.0, 5.0)

solucao = resolverSistema(penduloPID!, 
                          parametros_fisicos.estadosIniciais, 
                          intervaloTempo, 
                          parametros_controlador_PID)


plotarNoTempo(solucao, titulo="Estabilização PID Pêndulo", estados=1:3)
gerarAnimacao(solucao, parametros_fisicos, 30)