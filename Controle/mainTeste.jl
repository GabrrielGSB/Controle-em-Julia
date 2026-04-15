include("FuncSimuEDO.jl")
include("SistemasLTI.jl")

# paramSistema = PenduloParams(gravidade   = 9.81, 
#                              comprimento = 1.0, 
#                              massa       = 1.0, 
#                              coefAtrito  = 0.5,
#                              estadosIniciais = (θ = 2.0, ω = 0.0)) 

# paramSistema = MassaMolaParams(coefElastico    = 40.0, 
#                                coefAtrito      = 0.8, 
#                                massa           = 1.0, 
#                                estadosIniciais = (x = 0.5, v = 0.0))   

paramSistema = DroneParams(massa           = 0.468,
                           Ixx             = 4.856e-3,
                           Iyy             = 4.856e-3,
                           Izz             = 8.801e-3,
                           Ct              = 2.980e-6,
                           Cl              = 1.14e-7,
                           L               = 0.225,
                           controleW       = [800.0, 800.0, 810.0, 800.0],
                           estadosIniciais = (X=0.0,Y=0.0,Z=0.0,
                                              u=0.0,v=0.0,w=0.0,
                                              ϕ=0.5,θ=0.0,ψ=0.0,
                                              p=0.0,q=0.0,r=0.0),
                           )

simularSistemaMalhaAberta(sistema         = drone!, 
                          estadosIniciais = paramSistema.estadosIniciais, 
                          parametros      = paramSistema, 
                          tempoSimulacao  = (0.0, 10.0),
                          eixosVisu       = (0, 7),
                          visualizacao    = true,
                          visualizarTodosEstados = false,
                          animar          = false,
                          animarFPS       = 24)
