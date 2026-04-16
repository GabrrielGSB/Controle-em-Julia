using ProgressMeter
using Plots
include("../Métodos Controle/PID.jl")

#=
2. SISTEMA DINÂMICO DE PÊNDULO SIMPLES-----------------------------------------------------------------
    2.1. DESCRIÇÃO:
        Consiste em uma massa 'm' suspensa por uma haste rígida (ou fio) de 
        comprimento 'L', sob a influência da gravidade 'g'. O movimento é 
        oscilatório, onde a energia potencial gravitacional se converte em 
        energia cinética. Diferente do massa-mola, este sistema apresenta 
        uma dinâmica não-linear devido à componente restauradora da gravidade.

    2.2. MODELAGEM EM ESPAÇO DE ESTADOS:
        Baseado na 2ª Lei de Newton para rotações (τ = I*α), o sistema de 2ª 
        ordem é decomposto em duas equações de 1ª ordem:

        2.2.1. Vetor de Estados: 
            - x[1] : Ângulo θ (rad)        -> Posição angular
            - x[2] : Velocidade angular ω (rad/s)

        2.2.2. Vetor de Derivadas (dx):
            - dx[1] = x[2]          -> A variação do ângulo é a velocidade angular.
            - dx[2] = dω/dt         -> A variação da velocidade é a aceleração angular,
                                       calculada pelo torque gravitacional:
                                       α = -(g/L) * sen(x[1]) - (b/(m*L^2)) * x[2]

    2.3. PARÂMETROS:
        - g: Aceleração da gravidade [m/s²] (padrão: 9.81)
        - L: Comprimento da haste [m]
        - m: Massa da partícula [kg]
        - b: Coeficiente de atrito rotacional [N.m.s/rad] (opcional)
=#

struct PenduloParams
    gravidade       ::Float64
    comprimento     ::Float64
    massa           ::Float64
    coefAtrito      ::Float64
    estadosIniciais ::Vector{Float64}
    variaveisEstado ::Tuple{String, String}
    nomeSistema     ::String
    
    function PenduloParams(; gravidade, comprimento, massa, coefAtrito, estadosIniciais)
        new(gravidade, comprimento, massa, coefAtrito, estadosIniciais, 
            ("Ângulo θ (rad)", "Velocidade angular ω (rad/s)"), 
            "Pêndulo Simples")
    end
end

function pendulo!(dx, x, p, t)
    g, L, m, b = p.gravidade, p.comprimento, p.massa, p.coefAtrito

    x1, x2 = x[1], x[2]

    dx[1] = x2  
    dx[2] = -(g/L)*sin(x1) - (b/(m*L^2))*x2 
end
function penduloPID!(dx, x, p::PIDparams, t)
    g, L = p.fisica.gravidade, p.fisica.comprimento
    m, b = p.fisica.massa, p.fisica.coefAtrito
    
    I = m*L^2

    x1, x2 = x[1], x[2]

    u = calcularPID(x, p)

    dx[1] = x2  
    dx[2] = (-m*g*L*sin(x1) - b*x2 + u)/I
    dx[3] = p.referencia - x1
end

function calcularPID(x, p::PIDparams)
    θ             = x[1]
    ω             = x[2]
    integral_erro = x[3] 
    
    erro          = p.referencia - θ
    derivada_erro = -ω 

    return (p.Kp * erro) + (p.Ki * integral_erro) + (p.Kd * derivada_erro)
end

function gerarAnimacao(solucao, p::PenduloParams, fps)
    gr() 
    println("Iniciando animação do Pêndulo...")

    duracao = solucao.t[end]
    instantes_de_tempo = 0 : (1/fps) : duracao
    
    # 1. Criamos o objeto de animação vazio
    anim = Plots.Animation()
    
    println("Iniciando renderização do Pêndulo...")

    @showprogress "Progresso: " for t in instantes_de_tempo
        # Ângulo atual (interpolação contínua)
        θ = solucao(t)[1]
        L = p.comprimento
        
        # Trigonometria para posição da massa
        x_massa =  L * sin(θ)
        y_massa = -L * cos(θ)

        plot = Plots.plot([0, x_massa], [0, y_massa], 
                          lw=3, color=:black, label="", 
                          xlim=(-2, 2), ylim=(-1.2, 2),
                          aspect_ratio=:equal, title="Tempo: $(round(t, digits=2))s")
        
        scatter!(plot, [x_massa], [y_massa], 
                 markersize=10, color=:red, label="")
        
        Plots.frame(anim, plot)
    end
    
    # Configuração do caminho e salvamento
    caminho_video = "Controle/Animações/pendulo.mp4"
    mkpath(dirname(caminho_video))
    
    println("\nFinalizando arquivo de vídeo...")
    mp4(anim, caminho_video, fps = fps)
    println("Sucesso! Vídeo salvo em: $caminho_video")
end
