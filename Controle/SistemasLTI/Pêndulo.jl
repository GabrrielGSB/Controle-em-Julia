using ProgressMeter
import Plots

include("../Abstrações/Interfaces.jl")

#=
# PÊNDULO SIMPLES-----------------------------------------------------------------
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

# =========================================================================
# STRUCT DA PLANTA
    struct Pendulo <: SistemaPlanta
        comprimento     ::Float64
        massa           ::Float64
        coefAtrito      ::Float64
        estadosIniciais ::Vector{Float64}
        dim             ::Int
        dinamica!       ::Function               
        variaveisEstado ::Tuple{String, String}
        nomeSistema     ::String
        
        function Pendulo(; comprimento, massa, coefAtrito, estadosIniciais)
            new(comprimento, massa, coefAtrito, estadosIniciais,
                2,
                pendulo!,                        # ← aponta para a função abaixo
                ("Ângulo θ (rad)", "Velocidade angular ω (rad/s)"),
                "Pêndulo Simples")
        end
    end
# =========================================================================

# =========================================================================
# DINÂMICA DA PLANTA 
    """
        Dinâmica pura do pêndulo. O controle `u` é uma entrada externa opcional —
        por padrão zero (malha aberta).
    """
    function pendulo!(dx, x, p::Pendulo, t; u=0.0)   
        g = 9.81
        L, m, b = p.comprimento, p.massa, p.coefAtrito
        I = m * L^2

        dx[1] = x[2]
        dx[2] = (-m*g*L*sin(x[1]) - b*x[2] + u) / I
    end
# =========================================================================

# =========================================================================
# ANIMAÇÃO 
    function gerarAnimacao(solucao_raw, p::Pendulo, fps;
                           mostrar_controle::Bool = false)
        Plots.gr()
        println("Iniciando animação do Pêndulo...")

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

        seta_max           = 0.5
        duracao            = solucao.t[end]
        instantes_de_tempo = 0 : (1/fps) : duracao
        anim               = Plots.Animation()

        println("Iniciando renderização do Pêndulo...")

        @showprogress "Progresso: " for t in instantes_de_tempo
            θ = solucao(t)[1]
            ω = solucao(t)[2]
            L = p.comprimento

            x_massa =  L * sin(θ)
            y_massa = -L * cos(θ)

            plt = Plots.plot([0, x_massa], [0, y_massa],
                            lw=3, color=:black, label="",
                            xlim=(-2, 2), ylim=(-1.5, 1.5),
                            aspect_ratio=:equal,
                            title="Tempo: $(round(t, digits=2))s",
                            legend=:topright)

            Plots.scatter!(plt, [x_massa], [y_massa],
                        markersize=10, color=:red, label="Massa")

            if mostrar_controle
                idx = clamp(searchsortedlast(t_u, t), 1, length(u_historico))
                u   = u_historico[idx][1]

                tang_x  =  cos(θ)
                tang_y  =  sin(θ)
                escala  = (abs(u) / (u_max + 1e-10)) * seta_max
                seta_dx = sign(u) * tang_x * escala
                seta_dy = sign(u) * tang_y * escala

                Plots.quiver!(plt,
                            [x_massa], [y_massa],
                            quiver=([seta_dx], [seta_dy]),
                            color=:blue, lw=2.5, label="u(t)")

                info = "u  = $(round(u,          digits=2)) N·m\n" *
                    "θ  = $(round(rad2deg(θ), digits=1)) °\n"   *
                    "ω  = $(round(ω,          digits=3)) rad/s"

                Plots.annotate!(plt, -1.9, 1.4,
                                Plots.text(info, :left, 8, :darkblue))
            end

            Plots.frame(anim, plt)
        end

        caminho_video = "Controle/Animações/pendulo.mp4"
        mkpath(dirname(caminho_video))
        println("\nFinalizando arquivo de vídeo...")
        Plots.mp4(anim, caminho_video, fps=fps)
        println("Sucesso! Vídeo salvo em: $caminho_video")
    end
# =========================================================================