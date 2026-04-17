include("../../SistemasLTI/Drone.jl")
include("../../Funções Auxiliares/ResoluçãoEDO.jl")
include("../../Funções Auxiliares/Visualização.jl")


# Configurando os parâmetros
paramSistema = DroneParams(
    massa = 0.468,
    Ixx = 4.856e-3,
    Iyy = 4.856e-3,
    Izz = 8.801e-3,
    Ct = 2.980e-6,
    Cl = 1.14e-7,
    L = 0.225,
    estadosIniciais = zeros(12) # Drone começa no chão (Z=0) e parado
)

# Definindo o problema com a nova função
tempoSimulacao = (0.0, 25.0) 

solucao = resolverSistema(drone_malha_fechada!, 
                          paramSistema.estadosIniciais,
                          tempoSimulacao,
                          paramSistema)

plotarNoTempo(solucao, estados=(3,7,8))

# gerarAnimacao(solucao, paramSistema, 24)