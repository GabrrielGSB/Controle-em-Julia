#=
1. SISTEMA DINÂMICO MASSA-MOLA-AMORTECEDOR------------------------------------------------------------
    1.1. DESCRIÇÃO:
        Representa um corpo de massa 'm' preso a uma mola de rigidez elástica 'k' e um 
        amortecedor de coeficiente 'b'. O sistema converte energia potencial 
        (mola) em cinética (movimento) enquanto o amortecedor dissipa energia 
        na forma de calor (atrito).

    1.2. MODELAGEM EM ESPAÇO DE ESTADOS:
        O sistema de 2ª ordem (Lei de Newton: F = m*a) é decomposto em duas 
        equações de 1ª ordem para ser resolvido numericamente:

        1.2.1. Vetor de Variáveis de Estado: 
            - x[1] : Posição (m)
            - x[2] : Velocidade (m/s)

        1.2.2. Vetor de Derivadas (dx):
            - dx[1] = x[2]               -> A variação da posição é a velocidade.
            - dx[2] = dv/dt              -> A variação da velocidade é a aceleração, 

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