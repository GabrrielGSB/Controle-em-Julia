using Plots
using ProgressMeter

function gerarAnimacao(solucao, params, fps)
    println("Aviso: Nenhuma função de animação definida para o sistema: $(params.nomeSistema)")
end

#=
1. SISTEMA DINÂMICO MASSA-MOLA-AMORTECEDOR------------------------------------------------------------
    1.1. DESCRIÇÃO:
        Representa um corpo de massa 'm' preso a uma mola de rigidez 'k' e um 
        amortecedor de coeficiente 'b'. O sistema converte energia potencial 
        (mola) em cinética (movimento) enquanto o amortecedor dissipa energia 
        na forma de calor (atrito).

    1.2. MODELAGEM EM ESPAÇO DE ESTADOS:
        O sistema de 2ª ordem (Lei de Newton: F = m*a) é decomposto em duas 
        equações de 1ª ordem para ser resolvido numericamente:

        1.2.1. Vetor de Estados: 
            - x[1] : Posição (m)
            - x[2] : Velocidade (m/s)

        1.2.2. Vetor de Derivadas (dx):
            - dx[1] = x[2]         -> A variação da posição é a velocidade.
            - dx[2] = dv/dt        -> A variação da velocidade é a aceleração, 
                                        calculada pelo balanço de forças: 
                                        a = (-b*v - k*x) / m

    1.3.PARÂMETROS:
        - k: Constante elástica [N/m]
        - b: Coeficiente de amortecimento [N.s/m]
        - m: Massa do corpo [kg]
=#

struct MassaMolaParams
    coefElastico    ::Float64
    coefAtrito      ::Float64
    massa           ::Float64
    estadosIniciais ::NamedTuple
    variaveisEstado ::Tuple{String, String}
    nomeSistema     ::String

    function MassaMolaParams(; coefElastico, coefAtrito, massa, estadosIniciais)
        new(coefElastico, coefAtrito, massa, estadosIniciais, 
            ("Posição x (m)", "Velocidade v (m/s)"), 
            "Massa-Mola-Amortecedor")
    end
end

function massa_mola_amortecedor!(dx, x, params, t)
    k, b, m = params.coefElastico, params.coefAtrito, params.massa

    dx[1] = x[2]  
    dx[2] = (-b * x[2] - k * x[1]) / m  
end

function gerarAnimacao(solucao, params::MassaMolaParams, fps)
    gr() 
    
    duracao = solucao.t[end]
    instantes_de_tempo = 0 : (1/fps) : duracao
    
    # 1. Criamos um objeto de animação vazio
    anim = Plots.Animation()
    
    comprimento_repouso = 1.0  
    distancia_parede_min = 0.1 

    # 2. O @showprogress agora funciona sem conflitos!
    @showprogress "Renderizando Frames: " for t in instantes_de_tempo
        x_sol = solucao(t)[1]
        x_absoluto = comprimento_repouso + x_sol
        x_visual = max(distancia_parede_min, x_absoluto)
        
        # Geramos o gráfico do frame atual
        p = Plots.plot([0, x_visual], [0, 0], 
            lw=4, color=:blue, label="Mola", 
            xlim=(-0.2, 2.5), ylim=(-0.5, 0.5),
            aspect_ratio=:equal, grid=true)
            
        plot!(p, [0, 0], [-0.3, 0.3], color=:black, lw=6, label="")
        plot!(p, [-0.2, 2.5], [-0.12, -0.12], color=:grey, lw=1, label="")

        scatter!(p, [x_visual], [0], 
            markershape=:square, markersize=25, 
            color=:grey, label="Massa")
            
        title!(p, "Tempo: $(round(t, digits=2))s")

        Plots.frame(anim, p)
    end
    
    path_out = "Controle/Animações/massa_mola_reta.mp4"
    mkpath(dirname(path_out))
    
    println("\nCodificando vídeo final...")
    mp4(anim, path_out, fps = fps)
    println("Vídeo salvo em: $path_out")
end
#------------------------------------------------------------------------------------------------------

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
    estadosIniciais ::NamedTuple
    variaveisEstado ::Tuple{String, String}
    nomeSistema     ::String
    

    function PenduloParams(; gravidade, comprimento, massa, coefAtrito, estadosIniciais)
        new(gravidade, comprimento, massa, coefAtrito, estadosIniciais, 
            ("Ângulo θ (rad)", "Velocidade angular ω (rad/s)"), 
            "Pêndulo Simples")
    end
end

function pendulo!(dx, x, params, t)
        g, L, m, b = params.gravidade, params.comprimento, params.massa, params.coefAtrito

        dx[1] = x[2]  
        dx[2] = -(g / L) * sin(x[1]) - (b / (m * L^2)) * x[2]  
end

function gerarAnimacao(solucao, params::PenduloParams, fps)
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
        L = params.comprimento
        
        # Trigonometria para posição da massa
        x_massa =  L * sin(θ)
        y_massa = -L * cos(θ)

        # 3. Geramos o gráfico do frame atual e guardamos na variável 'p'
        p = Plots.plot([0, x_massa], [0, y_massa], 
            lw=3, color=:black, label="", 
            xlim=(-1.2, 1.2), ylim=(-1.2, 0.2),
            aspect_ratio=:equal, title="Tempo: $(round(t, digits=2))s")
        
        # Adicionamos a massa ao gráfico 'p'
        scatter!(p, [x_massa], [y_massa], 
                 markersize=10, color=:red, label="")
        
        # 4. Capturamos o frame manualmente
        Plots.frame(anim, p)
    end
    
    # Configuração do caminho e salvamento
    caminho_video = "Controle/Animações/pendulo.mp4"
    mkpath(dirname(caminho_video))
    
    println("\nFinalizando arquivo de vídeo...")
    mp4(anim, caminho_video, fps = fps)
    println("Sucesso! Vídeo salvo em: $caminho_video")
end

# function animacaoPendulo(; solucao, params, fps=30)
#     gr() # Backend otimizado para velocidade
    
#     duracao = solucao.t[end]
#     instantes_de_tempo = 0 : (1/fps) : duracao
    
#     # 1. Criamos o objeto de animação vazio
#     anim = Plots.Animation()
    
#     println("Iniciando renderização do Pêndulo...")

#     # 2. O @showprogress agora monitora o loop sem erros de macro
#     @showprogress "Progresso: " for t in instantes_de_tempo
#         # Ângulo atual (interpolação contínua)
#         θ = solucao(t)[1]
#         L = params.comprimento
        
#         # Trigonometria para posição da massa
#         x_massa =  L * sin(θ)
#         y_massa = -L * cos(θ)

#         # 3. Geramos o gráfico do frame atual e guardamos na variável 'p'
#         p = Plots.plot([0, x_massa], [0, y_massa], 
#             lw=3, color=:black, label="", 
#             xlim=(-1.2, 1.2), ylim=(-1.2, 0.2),
#             aspect_ratio=:equal, title="Tempo: $(round(t, digits=2))s")
        
#         # Adicionamos a massa ao gráfico 'p'
#         scatter!(p, [x_massa], [y_massa], 
#                  markersize=10, color=:red, label="")
        
#         # 4. Capturamos o frame manualmente
#         Plots.frame(anim, p)
#     end
    
#     # Configuração do caminho e salvamento
#     caminho_video = "Controle/Animações/pendulo.mp4"
#     mkpath(dirname(caminho_video))
    
#     println("\nFinalizando arquivo de vídeo...")
#     mp4(anim, caminho_video, fps = fps)
#     println("Sucesso! Vídeo salvo em: $caminho_video")
# end
#-------------------------------------------------------------------------------------------------------
