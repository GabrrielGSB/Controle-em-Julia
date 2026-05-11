using ProgressMeter
import Plots

include("../Abstrações/Interfaces.jl")

#=
# PÊNDULO INVERTIDO NO CARRINHO (CART-POLE) --------------------------------------
    1. DESCRIÇÃO:
        Representa um carrinho de massa 'M' que se move no eixo X 
        sob a ação de uma força 'u', tentando equilibrar uma haste de 
        comprimento 'l' e massa 'm' na ponta. É um clássico sistema subatuado
        (1 entrada de controle para 2 graus de liberdade físicos).

    2. MODELAGEM EM ESPAÇO DE ESTADOS:
        O sistema possui dinâmica não-linear acoplada entre a translação do
        carrinho e a rotação do pêndulo. 

        2.1. Vetor de Estados: 
            - x[1] : Posição do carrinho (m)
            - x[2] : Velocidade do carrinho (m/s)
            - x[3] : Ângulo do pêndulo θ (rad) -> 0 é para CIMA (invertido)
            - x[4] : Velocidade angular do pêndulo ω (rad/s)

        2.2. Vetor de Derivadas (dx):
            - dx[1] = x[2]
            - dx[2] = dv_c/dt (Aceleração linear do carrinho)
            - dx[3] = x[4]
            - dx[4] = dω/dt   (Aceleração angular do pêndulo)

    3. PARÂMETROS:
        - M: Massa do carrinho [kg] (padrão: 1.0)
        - m: Massa do pêndulo [kg] (padrão: 0.1)
        - l: Comprimento da haste [m] (padrão: 0.5)
        - g: Aceleração da gravidade [m/s²] (padrão: 9.81)
        - b: Atrito viscoso do carrinho [N.s/m] (padrão: 0.1)
=#

# =========================================================================
# STRUCT DA PLANTA
    struct PenduloInvertido <: Planta
        M               ::Float64
        m               ::Float64
        l               ::Float64
        g               ::Float64
        b               ::Float64
        estadosIniciais ::Vector{Float64}
        numEstados      ::Int
        variaveisEstado ::Tuple{String, String, String, String}
        nomeSistema     ::String
        dinamica!       ::Function 

        function PenduloInvertido(; 
                M = 1.0, 
                m = 0.1, 
                l = 0.5, 
                g = 9.81, 
                b = 0.1, 
                estadosIniciais = [0.0, 0.0, 0.1, 0.0]
            )
            new(M, m, l, g, b, estadosIniciais,
                4,
                ("Posição Carro (m)", "Velocidade Carro (m/s)", "Ângulo (rad)", "Velocidade Angular (rad/s)"),
                "Pêndulo Invertido no Carrinho",
                pendulo_invertido!)
        end
    end
# =========================================================================

# =========================================================================
# DINÂMICA DA PLANTA 
    """
        Dinâmica não-linear do Cart-Pole. O controle `u` atua linearmente
        sobre o carrinho.
    """
    function pendulo_invertido!(dx, x, p::PenduloInvertido, t; u=0.0)
        # Extração de parâmetros para limpeza das equações
        M, m, l, g, b = p.M, p.m, p.l, p.g, p.b
        
        # Extração dos estados
        x_c = x[1]
        v_c = x[2]
        θ   = x[3]
        ω   = x[4]
        
        # Termos trigonométricos recorrentes
        senθ = sin(θ)
        cosθ = cos(θ)
        
        # Denominador comum (matriz de inércia efetiva)
        D = M + m * (senθ^2)
        
        # Dinâmica translacional do carrinho
        dx[1] = v_c
        dx[2] = (u - b*v_c + m*l*(ω^2)*senθ - m*g*senθ*cosθ) / D
        
        # Dinâmica rotacional do pêndulo
        dx[3] = ω
        dx[4] = ((M + m)*g*senθ - cosθ*(u - b*v_c + m*l*(ω^2)*senθ)) / (l * D)
    end
# =========================================================================

# =========================================================================
# ANIMAÇÃO 
    function gerarAnimacao(solucao_raw, p::PenduloInvertido, fps;
                           mostrar_controle::Bool = false)
        Plots.gr()
        println("Iniciando animação do Cart-Pole...")

        # Desempacota dependendo do tipo de retorno de resolverSistema
        if mostrar_controle
            @assert solucao_raw isa NamedTuple "Para mostrar_controle=true, passe o retorno completo de resolverSistema(..., salvar_controle=true)"
            solucao     = solucao_raw.solucao
            t_u         = solucao_raw.t_u
            u_historico = solucao_raw.u
            u_max       = maximum(abs(v[1]) for v in u_historico)
        else
            solucao = solucao_raw isa NamedTuple ? solucao_raw.solucao : solucao_raw
        end

        seta_max           = 1.0
        duracao            = solucao.t[end]
        instantes_de_tempo = 0 : (1/fps) : duracao
        anim               = Plots.Animation()

        println("Iniciando renderização do Cart-Pole...")

        @showprogress "Progresso: " for t in instantes_de_tempo
            x_c = solucao(t)[1]
            θ   = solucao(t)[3]
            ω   = solucao(t)[4]
            l   = p.l

            # Coordenadas da massa do pêndulo (0 radianos = para cima)
            x_massa = x_c + l * sin(θ)
            y_massa = l * cos(θ)

            # Dimensões do carrinho para o desenho
            larg_carrinho = 0.6
            alt_carrinho  = 0.3
            cx = [x_c - larg_carrinho/2, x_c + larg_carrinho/2, x_c + larg_carrinho/2, x_c - larg_carrinho/2, x_c - larg_carrinho/2]
            cy = [-alt_carrinho/2, -alt_carrinho/2, alt_carrinho/2, alt_carrinho/2, -alt_carrinho/2]

            # Inicia o plot
            plt = Plots.plot(cx, cy, fill=true, color=:blue, label="Carrinho",
                             xlim=(-3, 3), ylim=(-1.5, 1.5),
                             aspect_ratio=:equal,
                             title="Tempo: $(round(t, digits=2))s",
                             legend=:topright)

            # Desenha o trilho
            Plots.plot!(plt, [-10, 10], [-alt_carrinho/2, -alt_carrinho/2], color=:black, lw=1, label="")

            # Desenha a haste
            Plots.plot!(plt, [x_c, x_massa], [0, y_massa], lw=3, color=:black, label="")

            # Desenha a massa do pêndulo
            Plots.scatter!(plt, [x_massa], [y_massa], markersize=8, color=:red, label="Massa")

            if mostrar_controle
                idx = clamp(searchsortedlast(t_u, t), 1, length(u_historico))
                u   = u_historico[idx][1]

                escala  = (abs(u) / (u_max + 1e-10)) * seta_max
                seta_dx = sign(u) * escala
                
                # Seta empurrando o carrinho horizontalmente
                Plots.quiver!(plt,
                              [x_c - sign(u)*larg_carrinho], [0.0],
                              quiver=([seta_dx], [0.0]),
                              color=:magenta, lw=3.0, label="u(t)")

                info = "u  = $(round(u,          digits=2)) N\n" *
                       "x  = $(round(x_c,        digits=2)) m\n" *
                       "θ  = $(round(rad2deg(θ), digits=1)) °"

                Plots.annotate!(plt, -2.5, 1.2,
                                Plots.text(info, :left, 8, :darkblue))
            end

            Plots.frame(anim, plt)
        end

        caminho_video = "Controle/Animações/pendulo_invertido.mp4"
        mkpath(dirname(caminho_video))
        println("\nFinalizando arquivo de vídeo...")
        Plots.mp4(anim, caminho_video, fps=fps)
        println("Sucesso! Vídeo salvo em: $caminho_video")
    end
# =========================================================================