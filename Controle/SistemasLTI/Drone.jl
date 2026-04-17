using ProgressMeter
using DynamicPolynomials

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
    estadosIniciais ::Vector{Float64}
    variaveisEstado ::Tuple
    nomeSistema     ::String

    function DroneParams(; massa, gravidade=9.81, Ixx, Iyy, Izz, Ct, Cl, L, controleW=[0.0, 0.0, 0.0, 0.0], estadosIniciais)
        new(massa, gravidade, Ixx, Iyy, Izz, Ct, Cl, L, controleW, estadosIniciais, 
            ("X", "Y", "Z", "Vx", "Vy", "Vz", "Φ", "θ", "Ψ", "VΦ", "Vθ", "VΨ"), 
            "Drone Quadrotor")
    end
end

function drone!(dx, x, p, t)
    m, g          = p.massa, p.gravidade
    Ixx, Iyy, Izz = p.Ixx, p.Iyy, p.Izz
    Ct, Cl, L     = p.Ct, p.Cl, p.L
    
    w1, w2, w3, w4 = p.controleW[1], p.controleW[2], p.controleW[3], p.controleW[4]

    # Cálculos das entradas de controle (U1 a U4) com base na velocidade dos motores
    U1 = Ct*(w1^2 + w2^2 + w3^2 + w4^2)
    U2 = Ct*L*(w1^2 - w3^2)
    U3 = Ct*L*(w2^2 - w4^2)
    U4 = Cl*(w1^2 - w2^2 + w3^2 - w4^2)

    # Extração de estados
    # X, Y, Z = x[1], x[2], x[3] 
    Vx, Vy, Vz = x[4], x[5], x[6] 
    Φ, θ, Ψ    = x[7], x[8], x[9]
    VΦ, Vθ, VΨ = x[10], x[11], x[12]

    # 1. Cinemática de Translação (Derivadas de X, Y e Z)
    dx[1] = Vz*(sin(Φ)*sin(Ψ) + cos(Φ)*cos(Ψ)*sin(θ)) - Vy*(cos(Φ)*sin(Ψ) - cos(Ψ)*sin(Φ)*sin(θ)) + Vx*cos(Ψ)*cos(θ)
    dx[2] = Vy*(cos(Φ)*cos(Ψ) + sin(Φ)*sin(Ψ)*sin(θ)) - Vz*(cos(Ψ)*sin(Φ) - cos(Φ)*sin(Ψ)*sin(θ)) + Vx*cos(θ)*sin(Ψ)
    dx[3] = Vz*cos(Φ)*cos(θ) - Vx*sin(θ) + Vy*cos(θ)*sin(Φ)

    # 2. Dinâmica de Translação (Derivadas de Vx, Vy e Vz)
    dx[4] = VΨ*Vy - Vθ*Vz + g*sin(θ)
    dx[5] = VΦ*Vz - VΨ*Vx - g*cos(θ)*sin(Φ)
    dx[6] = Vθ*Vx - VΦ*Vy - g*cos(Φ)*cos(θ) + (U1/m)

    # 3. Cinemática de Rotação (Derivadas de ϕ, θ e ψ)
    dx[7] = VΦ + Vθ*sin(Φ)*tan(θ) + VΨ*cos(Φ)*tan(θ)
    dx[8] = Vθ*cos(Φ) - VΨ*sin(Φ)
    dx[9] = VΨ*(cos(Φ)/cos(θ)) + Vθ*(sin(Φ)/cos(θ))

    # 4. Dinâmica de Rotação (VΦ, Vθ, VΨ dot)
    dx[10] = (1/Ixx) * (U2 + (Iyy - Izz)*Vθ*VΨ)
    dx[11] = (1/Iyy) * (U3 + (Izz - Ixx)*VΦ*VΨ)
    dx[12] = (1/Izz) * (U4 + (Ixx - Iyy)*VΦ*Vθ)
end

function drone_malha_fechada!(dx, x, p, t)
    m, g          = p.massa, p.gravidade
    Ixx, Iyy, Izz = p.Ixx, p.Iyy, p.Izz
    Ct, Cl, L     = p.Ct, p.Cl, p.L
    
    Vx, Vy, Vz = x[4], x[5], x[6] 
    Φ, θ, Ψ    = x[7], x[8], x[9]
    VΦ, Vθ, VΨ = x[10], x[11], x[12]
    Z = x[3] 

    # =========================================================
    # APLICAÇÃO DO CONTROLADOR
        # REFERÊNCIAS
            Z_ref = 1 
            Φ_ref = 0.087
            θ_ref = 0.087
            Ψ_ref = 0.0
        
        # CÁLCULO DOS ERROS
            erro_Z = Z - Z_ref
            erro_Φ = Φ - Φ_ref
            erro_θ = θ - θ_ref
            erro_Ψ = Ψ - Ψ_ref

        # DEFINIÇÃO DAS AÇÕES DE CONTROLE
            U_barra = -3 * erro_Z - 0.5 * Vz

            U1 = ((m*g) + U_barra) 
            U2 = -0.05*erro_Φ - 0.05*VΦ
            U3 = -0.05*erro_θ - 0.05*Vθ
            U4 = -0.0061*erro_Ψ - 0.002*VΨ
    # =========================================================

    # 1. Cinemática de Translação
    dx[1] = Vz*(sin(Φ)*sin(Ψ) + cos(Φ)*cos(Ψ)*sin(θ)) - Vy*(cos(Φ)*sin(Ψ) - cos(Ψ)*sin(Φ)*sin(θ)) + Vx*cos(Ψ)*cos(θ)
    dx[2] = Vy*(cos(Φ)*cos(Ψ) + sin(Φ)*sin(Ψ)*sin(θ)) - Vz*(cos(Ψ)*sin(Φ) - cos(Φ)*sin(Ψ)*sin(θ)) + Vx*cos(θ)*sin(Ψ)
    dx[3] = Vz*cos(Φ)*cos(θ) - Vx*sin(θ) + Vy*cos(θ)*sin(Φ)

    # 2. Dinâmica de Translação
    dx[4] = VΨ*Vy - Vθ*Vz + g*sin(θ)
    dx[5] = VΦ*Vz - VΨ*Vx - g*cos(θ)*sin(Φ)
    dx[6] = Vθ*Vx - VΦ*Vy - g*cos(Φ)*cos(θ) + (U1/m)

    # 3. Cinemática de Rotação
    dx[7] = VΦ + Vθ*sin(Φ)*tan(θ) + VΨ*cos(Φ)*tan(θ)
    dx[8] = Vθ*cos(Φ) - VΨ*sin(Φ)
    dx[9] = VΨ*(cos(Φ)/cos(θ)) + Vθ*(sin(Φ)/cos(θ))

    # 4. Dinâmica de Rotação
    dx[10] = (1/Ixx) * (U2 + (Iyy - Izz)*Vθ*VΨ)
    dx[11] = (1/Iyy) * (U3 + (Izz - Ixx)*VΦ*VΨ)
    dx[12] = (1/Izz) * (U4 + (Ixx - Iyy)*VΦ*Vθ)
end

function obterEqPolyDrone()
    @polyvar X Y Z u v w_lin phi theta psi p q r
    vars = [X, Y, Z, u, v, w_lin, phi, theta, psi, p, q, r]

    m, grav = 0.468, 9.81
    Ixx, Iyy, Izz = 4.856e-3, 4.856e-3, 8.801e-3

    s_phi = phi;     c_phi = 1 - 0.5*phi^2
    s_theta = theta; c_theta = 1 - 0.5*theta^2
    s_psi = psi;     c_psi = 1 - 0.5*psi^2
    sec_theta = 1 + 0.5*theta^2

    f = Vector{Any}(undef, 12)
    
    f[1] = w_lin * (s_phi * s_psi + c_phi * c_psi * s_theta) - v * (c_phi * s_psi - c_psi * s_phi * s_theta) + u * c_psi * c_theta
    f[2] = v * (c_phi * c_psi + s_phi * s_psi * s_theta) - w_lin * (c_psi * s_phi - c_phi * s_psi * s_theta) + u * c_theta * s_psi
    f[3] = w_lin * c_phi * c_theta - u * s_theta + v * c_theta * s_phi
    
    f[4] = r * v - q * w_lin + grav * s_theta
    f[5] = p * w_lin - r * u - grav * c_theta * s_phi
    f[6] = q * u - p * v + grav * (1 - c_phi * c_theta)
    
    f[7] = p + q * s_phi * theta + r * c_phi * theta
    f[8] = q * c_phi - r * s_phi
    f[9] = r * (c_phi * sec_theta) + q * (s_phi * sec_theta) 
    
    f[10] = ((Iyy - Izz) / Ixx) * q * r
    f[11] = ((Izz - Ixx) / Iyy) * p * r
    f[12] = ((Ixx - Iyy) / Izz) * p * q

    g = zeros(12, 4)
    g[6, 1]  = 1.0 / m
    g[10, 2] = 1.0 / Ixx
    g[11, 3] = 1.0 / Iyy
    g[12, 4] = 1.0 / Izz

    return f, g, vars
end


function gerarAnimacao(solucao, p::DroneParams, fps)
    gr() 
    println("Iniciando animação 3D do Drone...")

    duracao = solucao.t[end]
    instantes_de_tempo = 0 : (1/fps) : duracao
    
    anim = Plots.Animation()
    
    # Extrair limites baseados na trajetória total para fixar a câmera 3D
    traj_x = [solucao(t)[1] for t in solucao.t]
    traj_y = [solucao(t)[2] for t in solucao.t]
    traj_z = [solucao(t)[3] for t in solucao.t]
    
    # Adicionei uma margem um pouco maior para o drone não cortar nas bordas
    xlims = (minimum(traj_x)-1.5, maximum(traj_x)+1.5)
    ylims = (minimum(traj_y)-1.5, maximum(traj_y)+1.5)
    zlims = (minimum(traj_z)-1.5, maximum(traj_z)+1.5)

    @showprogress "Progresso: " for t in instantes_de_tempo
        # Extração de Posição
        X_atual = solucao(t)[1]
        Y_atual = solucao(t)[2]
        Z_atual = solucao(t)[3]
        
        # Extração de Atitude (Ângulos de Euler)
        Phi_atual   = solucao(t)[7]
        Theta_atual = solucao(t)[8]
        Psi_atual   = solucao(t)[9]
        
        # Histórico da trajetória até o tempo t
        hist_x = [solucao(tau)[1] for tau in 0:(1/fps):t]
        hist_y = [solucao(tau)[2] for tau in 0:(1/fps):t]
        hist_z = [solucao(tau)[3] for tau in 0:(1/fps):t]

        # 1. Desenhar o rastro da trajetória
        plt = Plots.plot(hist_x, hist_y, hist_z, 
            lw=2, color=:blue, label="Trajetória",
            xlims=xlims, ylims=ylims, zlims=zlims,
            xlabel="Eixo X", ylabel="Eixo Y", zlabel="Eixo Z",
            camera=(45, 30), 
            title="Posição 3D do Drone - $(round(t, digits=2))s")
        
        # =====================================================================
        # ADIÇÃO: CÁLCULO E DESENHO DA ESTRUTURA DO DRONE
        # =====================================================================
        
        # Tamanho visual do braço (ajuste o multiplicador se ficar muito pequeno)
        L_vis = p.L * 3.0 
        
        # Matrizes de Rotação (Z-Y-X)
        Rx = [1.0 0.0 0.0; 0.0 cos(Phi_atual) -sin(Phi_atual); 0.0 sin(Phi_atual) cos(Phi_atual)]
        Ry = [cos(Theta_atual) 0.0 sin(Theta_atual); 0.0 1.0 0.0; -sin(Theta_atual) 0.0 cos(Theta_atual)]
        Rz = [cos(Psi_atual) -sin(Psi_atual) 0.0; sin(Psi_atual) cos(Psi_atual) 0.0; 0.0 0.0 1.0]
        R = Rz * Ry * Rx # Rotação final combinada
        
        # Posições locais dos 4 motores (Configuração em X)
        m1_local = [ L_vis,  L_vis, 0.0] # Frente-Direita
        m2_local = [-L_vis, -L_vis, 0.0] # Trás-Esquerda
        m3_local = [ L_vis, -L_vis, 0.0] # Frente-Esquerda
        m4_local = [-L_vis,  L_vis, 0.0] # Trás-Direita
        
        # Rotacionando e transladando para o mundo real
        pos_centro = [X_atual, Y_atual, Z_atual]
        M1 = pos_centro + R * m1_local
        M2 = pos_centro + R * m2_local
        M3 = pos_centro + R * m3_local
        M4 = pos_centro + R * m4_local
        
        # Desenhando os braços do drone (Linhas pretas cruzadas)
        plot!(plt, [M2[1], M1[1]], [M2[2], M1[2]], [M2[3], M1[3]], color=:black, lw=3, label="")
        plot!(plt, [M4[1], M3[1]], [M4[2], M3[2]], [M4[3], M3[3]], color=:black, lw=3, label="")
        
        # Desenhando os 4 motores nas pontas (Bolinhas vermelhas)
        scatter!(plt, [M1[1], M2[1], M3[1], M4[1]], 
                      [M1[2], M2[2], M3[2], M4[2]], 
                      [M1[3], M2[3], M3[3], M4[3]], 
                      markershape=:circle, markersize=3, color=:red, label="")
        
        # Desenhando o centro de massa do drone (Bolinha verde central)
        scatter!(plt, [X_atual], [Y_atual], [Z_atual], 
                      markershape=:circle, markersize=5, color=:green, label="Drone")
        
        # =====================================================================

        Plots.frame(anim, plt)
    end
    
    caminho_video = "Controle/Animações/drone_3d.mp4"
    mkpath(dirname(caminho_video))
    
    println("\nFinalizando arquivo de vídeo 3D...")
    mp4(anim, caminho_video, fps = fps)
    println("Sucesso! Vídeo salvo em: $caminho_video")
end