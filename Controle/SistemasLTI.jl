using Plots
using ProgressMeter

function gerarAnimacao(solucao, params, fps)
    println("Aviso: Nenhuma função de animação definida para o sistema: $(params.nomeSistema)")
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

#=
3. SISTEMA DINÂMICO DE DRONE (QUADROTOR 6-DOF) --------------------------------------------------------
    3.1. DESCRIÇÃO:
        Representa a cinemática e a dinâmica não-linear de um drone no espaço 3D.
        O modelo considera forças e torques atuando no referencial do corpo e os 
        converte para o referencial de inércia (Terra).

    3.2. MODELAGEM EM ESPAÇO DE ESTADOS:
        O sistema possui 12 variáveis de estado, decompostas em 12 equações de 1ª ordem.

        3.2.1. Vetor de Estados: 
            - x[1], x[2], x[3]   : Posição na Terra (X, Y, Z)
            - x[4], x[5], x[6]   : Velocidade linear no eixo do corpo (u, v, w)
            - x[7], x[8], x[9]   : Ângulos de Euler (ϕ Roll, θ Pitch, ψ Yaw)
            - x[10], x[11], x[12]: Velocidade angular no eixo do corpo (p, q, r)

        3.2.2. Variáveis de Entrada de Controle (U):
            - U1 : Empuxo Total (Thrust)
            - U2 : Torque de Rolagem (Roll)
            - U3 : Torque de Arfagem (Pitch)
            - U4 : Torque de Guinada (Yaw)

    3.3. PARÂMETROS:
        - m: Massa do drone [kg]
        - g: Aceleração da gravidade [m/s²]
        - Ixx, Iyy, Izz: Momentos de inércia nos respectivos eixos [kg.m²]
        - controleU: Vetor com as entradas [U1, U2, U3, U4]
=#

struct DroneParams
    massa           ::Float64
    gravidade       ::Float64
    Ixx             ::Float64
    Iyy             ::Float64
    Izz             ::Float64
    Ct              ::Float64 
    Cl              ::Float64 
    L               ::Float64 
    controleW       ::Vector{Float64} # [w1, w2, w3, w4]
    estadosIniciais ::NamedTuple
    variaveisEstado ::Tuple
    nomeSistema     ::String

    function DroneParams(; massa, gravidade=9.81, Ixx, Iyy, Izz, Ct, Cl, L, controleW=[0.0, 0.0, 0.0, 0.0], estadosIniciais)
        new(massa, gravidade, Ixx, Iyy, Izz, Ct, Cl, L, controleW, estadosIniciais, 
            ("X", "Y", "Z", "u", "v", "w", "ϕ", "θ", "ψ", "p", "q", "r"), 
            "Drone Quadrotor")
    end
end

function drone!(dx, x, params, t)
    # Extração de parâmetros físicos
    m, g = params.massa, params.gravidade
    Ixx, Iyy, Izz = params.Ixx, params.Iyy, params.Izz
    Ct, Cl, L = params.Ct, params.Cl, params.L
    
    # Extração das velocidades dos motores
    w1, w2, w3, w4 = params.controleW[1], params.controleW[2], params.controleW[3], params.controleW[4]

    # Cálculos das entradas de controle (U1 a U4) com base na velocidade dos motores
    U1 = Ct * (w1^2 + w2^2 + w3^2 + w4^2)
    U2 = Ct * L * (w1^2 - w3^2)
    U3 = Ct * L * (w2^2 - w4^2)
    U4 = Cl * (w1^2 - w2^2 + w3^2 - w4^2)

    # Extração de estados
    # X, Y, Z = x[1], x[2], x[3] (Não são diretamente necessários para calcular as derivadas)
    u, v, w_lin     = x[4], x[5], x[6] # w_lin para não conflitar com velocidade dos motores
    phi, theta, psi = x[7], x[8], x[9]
    p, q, r         = x[10], x[11], x[12]

    # Pré-cálculo trigonométrico
    s_phi, c_phi = sin(phi), cos(phi)
    s_theta, c_theta, t_theta = sin(theta), cos(theta), tan(theta)
    s_psi, c_psi = sin(psi), cos(psi)

    # 1. Cinemática de Translação (X, Y, Z dot)
    dx[1] = w_lin * (s_phi * s_psi + c_phi * c_psi * s_theta) - v * (c_phi * s_psi - c_psi * s_phi * s_theta) + u * c_psi * c_theta
    dx[2] = v * (c_phi * c_psi + s_phi * s_psi * s_theta) - w_lin * (c_psi * s_phi - c_phi * s_psi * s_theta) + u * c_theta * s_psi
    dx[3] = w_lin * c_phi * c_theta - u * s_theta + v * c_theta * s_phi

    # 2. Dinâmica de Translação (u, v, w dot)
    dx[4] = r * v - q * w_lin + g * s_theta
    dx[5] = p * w_lin - r * u - g * c_theta * s_phi
    dx[6] = q * u - p * v - g * c_phi * c_theta + (U1 / m)

    # 3. Cinemática de Rotação (ϕ, θ, ψ dot)
    dx[7] = p + q * s_phi * t_theta + r * c_phi * t_theta
    dx[8] = q * c_phi - r * s_phi
    dx[9] = r * (c_phi / c_theta) + q * (s_phi / c_theta)

    # 4. Dinâmica de Rotação (p, q, r dot)
    dx[10] = (1 / Ixx) * (U2 + (Iyy - Izz) * q * r)
    dx[11] = (1 / Iyy) * (U3 + (Izz - Ixx) * p * r)
    dx[12] = (1 / Izz) * (U4 + (Ixx - Iyy) * p * q)
end

function gerarAnimacao(solucao, params::DroneParams, fps)
    gr() 
    println("Iniciando animação 3D do Drone...")

    duracao = solucao.t[end]
    instantes_de_tempo = 0 : (1/fps) : duracao
    
    anim = Plots.Animation()
    
    # Extrair limites baseados na trajetória total para fixar a câmera 3D
    traj_x = [solucao(t)[1] for t in solucao.t]
    traj_y = [solucao(t)[2] for t in solucao.t]
    traj_z = [solucao(t)[3] for t in solucao.t]
    
    xlims = (minimum(traj_x)-1, maximum(traj_x)+1)
    ylims = (minimum(traj_y)-1, maximum(traj_y)+1)
    zlims = (minimum(traj_z)-1, maximum(traj_z)+1)

    @showprogress "Progresso 3D: " for t in instantes_de_tempo
        X_atual = solucao(t)[1]
        Y_atual = solucao(t)[2]
        Z_atual = solucao(t)[3]
        
        # Histórico da trajetória até o tempo t
        hist_x = [solucao(tau)[1] for tau in 0:(1/fps):t]
        hist_y = [solucao(tau)[2] for tau in 0:(1/fps):t]
        hist_z = [solucao(tau)[3] for tau in 0:(1/fps):t]

        # Desenhar o rastro
        p = Plots.plot(hist_x, hist_y, hist_z, 
            lw=2, color=:blue, label="Trajetória",
            xlims=xlims, ylims=ylims, zlims=zlims,
            xlabel="Eixo X", ylabel="Eixo Y", zlabel="Eixo Z",
            camera=(45, 30), # Ângulo isométrico de visualização
            title="Posição 3D do Drone - $(round(t, digits=2))s")
        
        # Desenhar a posição da massa do drone
        scatter!(p, [X_atual], [Y_atual], [Z_atual], 
            markershape=:circle, markersize=6, color=:red, label="Drone")

        Plots.frame(anim, p)
    end
    
    caminho_video = "Controle/Animações/drone_3d.mp4"
    mkpath(dirname(caminho_video))
    
    println("\nFinalizando arquivo de vídeo 3D...")
    mp4(anim, caminho_video, fps = fps)
    println("Sucesso! Vídeo salvo em: $caminho_video")
end