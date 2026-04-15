# include("../SistemasLTI/Pêndulo.jl")
include("../Funções Auxiliares/Visualização.jl")
include("../Funções Auxiliares/ResoluçãoEDO.jl")

struct PenduloParams
    gravidade       ::Float64
    comprimento     ::Float64
    massa           ::Float64
    coefAtrito      ::Float64
    estadosIniciais ::NamedTuple
    variaveisEstado ::Tuple{String, String}
    nomeSistema     ::String
    

    function PenduloParams(; gravidade, comprimento, massa, coefAtrito, estadosIniciais)
        new(gravidade, comprimento, massa, coefAtrito, estadosIniciais, 
            ("Ângulo θ (rad)", "Velocidade angular ω (rad/s)"), 
            "Pêndulo Simples")
    end
end

# Novo parâmetro que abraça a física e o controle
struct PenduloPIDParams
    fisica::PenduloParams
    Kp::Float64
    Ki::Float64
    Kd::Float64
    setpoint::Float64 # Ângulo alvo em radianos
end

# Lógica do PID pura (Fora da função do sistema)
function calcular_torque_pid(x, p::PenduloPIDParams)
    θ = x[1]
    ω = x[2]
    erro_integral = x[3] # O terceiro estado é a integral do erro
    
    erro = p.setpoint - θ
    derro = -ω # Derivada da referência(0) - derivada do estado(ω)
    
    # Lei de controle PID
    u = (p.Kp * erro) + (p.Ki * erro_integral) + (p.Kd * derro)
    return u
end

function sistema_controlado!(dx, x, p::PenduloPIDParams, t)
    # 1. Recuperar parâmetros físicos
    params = p.fisica
    g, L, m, b = params.gravidade, params.comprimento, params.massa, params.coefAtrito
    I = m * L^2 # Momento de inércia para massa pontual
    
    # 2. Calcular Torque de Controle (Chamando a lógica externa)
    u = calcular_torque_pid(x, p)
    
    # 3. Equações de Estado
    # dx[1]: Variação do ângulo (Velocidade)
    dx[1] = x[2] 
    
    # dx[2]: Variação da velocidade (Aceleração angular)
    # Torque Total = Torque Gravitacional + Torque Atrito + Torque Controle
    torque_grav = -m * g * L * sin(x[1])
    torque_atrito = -b * x[2]
    
    dx[2] = (torque_grav + torque_atrito + u) / I
    
    # dx[3]: Variação da integral do erro (O próprio erro atual)
    dx[3] = p.setpoint - x[1]
end

# 1. Definir a Física
fisica = PenduloParams(
    gravidade=9.81, 
    comprimento=1.0, 
    massa=1.0, 
    coefAtrito=0.1, 
    estadosIniciais=(θ=0.0, ω=0.0) # Agora é uma NamedTuple!
)

# 2. Configurar o PID (Ajuste esses valores!)
# Para o pêndulo invertido (setpoint = π), os ganhos precisam ser altos.
params_com_controle = PenduloPIDParams(
    fisica, 
    50.0,   # Kp
    10.0,   # Ki
    15.0,   # Kd
    3.1415  # Setpoint (Vertical para cima)
)

# 3. Condições iniciais (Ângulo, Velocidade, Erro Integral Inicial)
x0 = [2.8, 0.0, 0.0] # Começa quase no topo (2.8 rad)
tspan = (0.0, 10.0)

# 4. Resolver (Usando sua função de ResoluçãoEDO.jl)
sol = resolverSistema(sistema_controlado!, x0, tspan, params_com_controle)

# 5. Visualizar (Usando sua função de Visualização.jl)
# Plotamos os estados 1 (ângulo) e 2 (velocidade)
plotarNoTempo(sol, titulo="Estabilização PID Pêndulo", estados=1:2)