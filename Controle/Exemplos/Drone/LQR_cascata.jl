
# =========================================================
# INCLUDES
    using LinearAlgebra
    using PlotlyJS
    using Printf

    include("../../SistemasLTI/Drone.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

# =========================================================
# SOLVER: EQUAÇÃO ALGÉBRICA DE RICCATI CONTÍNUA (CARE)
    function resolverCARE(A::Matrix{Float64}, B::Matrix{Float64},
                        Q::Matrix{Float64}, R::Matrix{Float64})
        n    = size(A, 1)
        Rinv = inv(R)
        H = [ A           -B*Rinv*B';
            -Q           -A'       ]

        vals, vecs = eigen(H)
        ordem   = sortperm(real.(vals))
        estavel = ordem[1:n]

        if any(real.(vals[estavel]) .≥ -1e-9)
            error("resolverCARE: não foi possível encontrar autovalores estáveis.")
        end

        V  = vecs[:, estavel]
        V1 = V[1:n, :]
        V2 = V[n+1:end, :]

        P = real(V2 / V1)
        P = Matrix(Symmetric((P + P') / 2))

        K     = Rinv * B' * P
        polos = eigvals(A - B * K)

        return K, P, polos
    end
# =========================================================

# =========================================================
# STRUCTS DOS SUBSISTEMAS
    struct SubsistemaLQR
        A     ::Matrix{Float64}
        B     ::Matrix{Float64}
        Q     ::Matrix{Float64}
        R     ::Matrix{Float64}
        K     ::Matrix{Float64}
        P     ::Matrix{Float64}
        polos ::Vector{ComplexF64}
        nome  ::String
    end

    function SubsistemaLQR(A::Matrix{Float64}, B::Matrix{Float64},
                        Q::Matrix{Float64}, R::Matrix{Float64};
                        nome::String = "Subsistema")
        K, P, polos = resolverCARE(A, B, Q, R)
        return SubsistemaLQR(A, B, Q, R, K, P, polos, nome)
    end

    struct ControladorAtitudeLQR
        roll   ::SubsistemaLQR
        pitch  ::SubsistemaLQR
        yaw    ::SubsistemaLQR
        lim_U2 ::Float64
        lim_U3 ::Float64
        lim_U4 ::Float64
    end

    struct ControladorAltitudeLQR
        sub      ::SubsistemaLQR
        massa    ::Float64
        lim_min  ::Float64
        lim_max  ::Float64
    end

    struct ControladorPosicaoLQR
        sub_X      ::SubsistemaLQR
        sub_Y      ::SubsistemaLQR
        lim_angulo ::Float64
        gravidade  ::Float64
    end

    struct ControladorCascataLQR
        posicao  ::ControladorPosicaoLQR
        altitude ::ControladorAltitudeLQR
        atitude  ::ControladorAtitudeLQR
    end
# =========================================================

# =========================================================
# CONSTRUTORES DOS CONTROLADORES
    function construirAtitudeLQR(Ixx::Float64, Iyy::Float64, Izz::Float64;
                                Q_roll  = [400.0 0.0; 0.0  4.0], R_roll  = reshape([4.0], 1, 1),
                                Q_pitch = [400.0 0.0; 0.0  4.0], R_pitch = reshape([4.0], 1, 1),
                                Q_yaw   = [200.0 0.0; 0.0  2.0], R_yaw   = reshape([8.0], 1, 1),
                                lim_U2 = 0.60, lim_U3 = 0.60, lim_U4 = 0.25)
        A = [0.0  1.0; 0.0  0.0]
        B_roll  = reshape([0.0;  1.0/Ixx], 2, 1)
        B_pitch = reshape([0.0;  1.0/Iyy], 2, 1)
        B_yaw   = reshape([0.0;  1.0/Izz], 2, 1)

        return ControladorAtitudeLQR(
            SubsistemaLQR(A, B_roll,  Q_roll,  R_roll;  nome = "Roll  Φ → U2"),
            SubsistemaLQR(A, B_pitch, Q_pitch, R_pitch; nome = "Pitch θ → U3"),
            SubsistemaLQR(A, B_yaw,   Q_yaw,   R_yaw;   nome = "Yaw   Ψ → U4"),
            lim_U2, lim_U3, lim_U4
        )
    end

    function construirAltitudeLQR(massa::Float64; gravidade = 9.81,
                                Q_Z = [200.0  0.0; 0.0  50.0], R_Z = reshape([1.0], 1, 1),
                                lim_min = 0.0, lim_max = 2.5 * massa * 9.81)
        A = [0.0  1.0; 0.0  0.0]
        B = reshape([0.0;  1.0/massa], 2, 1)
        return ControladorAltitudeLQR(SubsistemaLQR(A, B, Q_Z, R_Z; nome = "Altitude Z → U1"), massa, lim_min, lim_max)
    end

    function construirPosicaoLQR(; gravidade = 9.81,
                                Q_X = [50.0  0.0; 0.0  20.0], R_X = reshape([1.0], 1, 1),
                                Q_Y = [50.0  0.0; 0.0  20.0], R_Y = reshape([1.0], 1, 1),
                                lim_angulo = 0.35)
        A = [0.0  1.0; 0.0  0.0]
        B = reshape([0.0;  1.0], 2, 1)
        return ControladorPosicaoLQR(
            SubsistemaLQR(A, B, Q_X, R_X; nome = "Posição X → θ_des"),
            SubsistemaLQR(A, B, Q_Y, R_Y; nome = "Posição Y → Φ_des"),
            lim_angulo, gravidade
        )
    end
# =========================================================

# =========================================================
# PARÂMETROS DA SIMULAÇÃO
    struct DroneTrajetoriaLQRParams
        drone      ::Drone
        ctrl       ::ControladorCascataLQR
        trajetoria ::Function
    end
# =========================================================

# =========================================================
# DINÂMICA DA MALHA FECHADA EM CASCATA 
    function drone_lqr_cascata_mf!(dx, x, p::DroneTrajetoriaLQRParams, t)
        d  = p.drone
        c  = p.ctrl
        g  = d.gravidade

        # ── Extração dos 12 estados do drone ─────────────────────────
        X,  Y,  Z  = x[1], x[2], x[3]
        Vx, Vy, Vz = x[4], x[5], x[6]
        Φ,  θ,  Ψ  = x[7], x[8], x[9]
        VΦ, Vθ, VΨ = x[10], x[11], x[12]

        # ── Referência de trajetória com Velocidades ──────────────────
        ref = p.trajetoria(t)
        X_r, Y_r, Z_r, Ψ_r = ref.X, ref.Y, ref.Z, ref.Ψ
        
        # Extrai as velocidades de referência (se a trajetória fornecer)
        Vx_r = hasproperty(ref, :Vx) ? ref.Vx : 0.0
        Vy_r = hasproperty(ref, :Vy) ? ref.Vy : 0.0
        Vz_r = hasproperty(ref, :Vz) ? ref.Vz : 0.0
        VΨ_r = hasproperty(ref, :VΨ) ? ref.VΨ : 0.0

        # ══════════════════════════════════════════════════════════════
        # NÍVEL 2 — ALTITUDE (Rastreamento de Erro LQR: Δx = x - x_r)
        # ══════════════════════════════════════════════════════════════
        e_alt = [Z - Z_r;  Vz - Vz_r]
        δU1   = -(c.altitude.sub.K * e_alt)[1]
        U1    = clamp(c.altitude.massa * g + δU1, c.altitude.lim_min, c.altitude.lim_max)

        # ══════════════════════════════════════════════════════════════
        # NÍVEL 1 — POSIÇÃO HORIZONTAL
        # ══════════════════════════════════════════════════════════════
        Vx_i =  Vx * cos(Ψ) - Vy * sin(Ψ)
        Vy_i =  Vx * sin(Ψ) + Vy * cos(Ψ)

        # LQR atuando sobre o erro absoluto de estado e velocidade
        e_pos_X = [X - X_r;  Vx_i - Vx_r]
        e_pos_Y = [Y - Y_r;  Vy_i - Vy_r]

        ax_des = -(c.posicao.sub_X.K * e_pos_X)[1]
        ay_des = -(c.posicao.sub_Y.K * e_pos_Y)[1]

        θ_des = clamp(( cos(Ψ) * ax_des + sin(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)
        Φ_des = clamp(( sin(Ψ) * ax_des - cos(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)

        # ══════════════════════════════════════════════════════════════
        # NÍVEL 3 — ATITUDE 
        # ══════════════════════════════════════════════════════════════
        # Referências de taxa de Roll/Pitch assumidas 0 (pois são atuadores para a posição)
        e_roll  = [Φ - Φ_des;  VΦ - 0.0]
        e_pitch = [θ - θ_des;  Vθ - 0.0]
        e_yaw   = [atan(sin(Ψ - Ψ_r), cos(Ψ - Ψ_r));  VΨ - VΨ_r]

        U2 = clamp(-(c.atitude.roll.K  * e_roll )[1], -c.atitude.lim_U2, c.atitude.lim_U2)
        U3 = clamp(-(c.atitude.pitch.K * e_pitch)[1], -c.atitude.lim_U3, c.atitude.lim_U3)
        U4 = clamp(-(c.atitude.yaw.K   * e_yaw  )[1], -c.atitude.lim_U4, c.atitude.lim_U4)

        d.dinamica!(dx, x, d, t; u = [U1, U2, U3, U4])
    end
# =========================================================

# =========================================================
# TRAJETÓRIAS DE REFERÊNCIA COM DERIVADA ANALÍTICA (VELOCIDADES)
    function trajetoriaHover(t; altitude=2.0, t_subida=5.0)
        τ = clamp(t / t_subida, 0.0, 1.0)
        h = altitude * τ^2 * (3 - 2τ)
        vh = (0.0 <= t < t_subida) ? altitude * (6τ - 6τ^2) / t_subida : 0.0
        return (X=0.0, Y=0.0, Z=h, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=vh, VΨ=0.0)
    end

    function trajetoriaHelice(t; raio=1.5, ω_giro=0.4, v_subida=0.10, t_espera=5.0)
        if t < t_espera
            τ = t / t_espera
            h = τ^2 * (3 - 2τ)
            vh = (0.0 <= t < t_espera) ? (6τ - 6τ^2) / t_espera : 0.0
            return (X=0.0, Y=0.0, Z=h, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=vh, VΨ=0.0)
        end
        τ = t - t_espera
        return (X  = raio * cos(ω_giro * τ),
                Y  = raio * sin(ω_giro * τ),
                Z  = 1.0 + v_subida * τ,
                Ψ  = 0.0,
                Vx = -raio * ω_giro * sin(ω_giro * τ),
                Vy =  raio * ω_giro * cos(ω_giro * τ),
                Vz = v_subida,
                VΨ = 0.0)
    end

    function trajetoriaFigura8(t; escala=2.0, T=20.0, altitude=2.0, t_subida=5.0)
        if t < t_subida
            τ = t / t_subida
            vh = (0.0 <= t < t_subida) ? altitude * (6τ - 6τ^2) / t_subida : 0.0
            return (X=0.0, Y=0.0, Z=altitude * τ^2 * (3-2τ), Ψ=0.0, Vx=0.0, Vy=0.0, Vz=vh, VΨ=0.0)
        end
        τ = t - t_subida;  ω = 2π/T
        X = escala * sin(ω * τ)
        Y = escala * sin(ω * τ) * cos(ω * τ)
        Vx = escala * ω * cos(ω * τ)
        Vy = escala * ω * (cos(ω * τ)^2 - sin(ω * τ)^2)
        return (X=X, Y=Y, Z=altitude, Ψ=0.0, Vx=Vx, Vy=Vy, Vz=0.0, VΨ=0.0)
    end

    function trajetoriaWaypoints(t, waypoints; T_seg=5.0)
        n   = size(waypoints, 1)
        idx = clamp(floor(Int, t / T_seg) + 1, 1, n-1)
        τ   = clamp((t - (idx-1)*T_seg) / T_seg, 0.0, 1.0)
        s   = τ^2 * (3 - 2τ)
        ds  = (0.0 <= τ <= 1.0) ? (6τ - 6τ^2) / T_seg : 0.0
        
        p0  = waypoints[idx, :];  p1 = waypoints[min(idx+1, n), :]
        r   = p0 + s * (p1 - p0)
        v   = ds * (p1 - p0)
        
        return (X=r[1], Y=r[2], Z=r[3], Ψ=r[4], Vx=v[1], Vy=v[2], Vz=v[3], VΨ=v[4])
    end
# =========================================================

# =========================================================
# VISUALIZAÇÃO
    function imprimirRelatorioLQR(ctrl::ControladorCascataLQR)
        subsistemas = [ctrl.atitude.roll, ctrl.atitude.pitch, ctrl.atitude.yaw,
                    ctrl.altitude.sub, ctrl.posicao.sub_X, ctrl.posicao.sub_Y]

        println("\n╔══════════════════════════════════════════════════════════════════╗")
        println("║            RELATÓRIO DOS CONTROLADORES LQR EM CASCATA           ║")
        println("╠══════════════════════════════════════════════════════════════════╣")

        for sub in subsistemas
            K1, K2 = sub.K[1, 1], sub.K[1, 2]
            p1, p2 = sub.polos[1], sub.polos[2]
            est = all(real.(sub.polos) .< 0) ? "✓ ESTÁVEL" : "✗ INSTÁVEL"
            println()
            @printf("  ┌─ %s\n", sub.nome)
            @printf("  │  Ganho  K = [ %8.4f  %8.4f ]\n", K1, K2)
            @printf("  │  Polo 1   = %+.4f %+.4f·j   |  ωₙ ≈ %.2f rad/s\n", real(p1), imag(p1), abs(p1))
            @printf("  │  Polo 2   = %+.4f %+.4f·j   |  ζ  ≈ %.3f\n", real(p2), imag(p2), -real(p1) / (abs(p1) + 1e-12))
            @printf("  └─ %s\n", est)
        end
        println("\n╚══════════════════════════════════════════════════════════════════╝\n")
    end

    function plotarTrajetoria3D(sol, trajetoria, tspan; titulo = "LQR Cascata — Trajetória 3D")
        t_vec = sol.t
        refs  = [trajetoria(t) for t in t_vec]
        X_r = [r.X for r in refs];  Y_r = [r.Y for r in refs];  Z_r = [r.Z for r in refs]
        X   = [u[1] for u in sol.u]; Y = [u[2] for u in sol.u]; Z = [u[3] for u in sol.u]

        display(plot([
            scatter3d(x=X_r, y=Y_r, z=Z_r, mode="lines", name="Referência", line=attr(color="#d62728", width=3, dash="dash")),
            scatter3d(x=X,   y=Y,   z=Z,   mode="lines", name="Drone (real)", line=attr(color="#1f77b4", width=4)),
        ], Layout(title_text=titulo, title_x=0.5, scene=attr(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)"), width=850, height=650)))
    end

    function plotarRespostaCompleta(sol, trajetoria; titulo="LQR Cascata — Resposta Completa")
        t    = sol.t
        refs = [trajetoria(τ) for τ in t]

        X_r = [r.X for r in refs];  Y_r = [r.Y for r in refs]
        Z_r = [r.Z for r in refs];  Ψ_r = [r.Ψ for r in refs]

        X   = [u[1]  for u in sol.u];  Y   = [u[2]  for u in sol.u]
        Z   = [u[3]  for u in sol.u];  Φ   = [u[7]  for u in sol.u]
        θ   = [u[8]  for u in sol.u];  Ψ   = [u[9]  for u in sol.u]

        fig = make_subplots(rows=3, cols=2, shared_xaxes=true, vertical_spacing=0.08, horizontal_spacing=0.10,
                            subplot_titles=reshape(["Altitude Z (m)", "Posição X (m)", "Roll Φ (°)", "Posição Y (m)", "Pitch θ (°)", "Yaw Ψ (°)"], 1, 6))

        for (linha, (dados, ref, cor, nome)) in enumerate([(Z, Z_r, "#1f77b4", "Z"), (rad2deg.(Φ), zeros(length(t)), "#d62728", "Φ"), (rad2deg.(θ), zeros(length(t)), "#ff7f0e", "θ")])
            add_trace!(fig, scatter(x=t, y=ref, mode="lines", name="Ref $nome", line=attr(color="black", width=1, dash="dot"), showlegend=(linha==1)), row=linha, col=1)
            add_trace!(fig, scatter(x=t, y=dados, mode="lines", name=nome, line=attr(color=cor, width=2)), row=linha, col=1)
        end
        for (linha, (dados, ref, cor, nome)) in enumerate([(X, X_r, "#2ca02c", "X"), (Y, Y_r, "#9467bd", "Y"), (rad2deg.(Ψ), rad2deg.(Ψ_r), "#8c564b", "Ψ")])
            add_trace!(fig, scatter(x=t, y=ref, mode="lines", name="Ref $nome", line=attr(color="black", width=1, dash="dot"), showlegend=(linha==1)), row=linha, col=2)
            add_trace!(fig, scatter(x=t, y=dados, mode="lines", name=nome, line=attr(color=cor, width=2)), row=linha, col=2)
        end
        relayout!(fig, title_text=titulo, title_x=0.5, height=750, width=1000, hovermode="x unified")
        display(fig)
    end
# =========================================================

# =========================================================
# CONFIGURAÇÃO E EXECUÇÃO DO DRONE
    drone = Drone(massa = 0.468, gravidade = 9.81, Ixx = 4.856e-3, 
                Iyy = 4.856e-3, Izz = 8.801e-3, Ct = 2.980e-6, 
                Cl = 1.14e-7, L = 0.225, estadosIniciais = zeros(12))

    # Ganhos LQR otimizados
    ctrl_atitude = construirAtitudeLQR(drone.Ixx, drone.Iyy, drone.Izz;
        Q_roll  = [1028.4 0.0; 0.0 639.3], R_roll  = reshape([889.30], 1, 1),
        Q_pitch = [1715.3 0.0; 0.0 18.2], R_pitch = reshape([290.02], 1, 1),
        Q_yaw   = [551.1 0.0; 0.0 1393.3], R_yaw   = reshape([569.36], 1, 1))

    ctrl_altitude = construirAltitudeLQR(drone.massa; gravidade=drone.gravidade,
        Q_Z = [787.9 0.0; 0.0 796.9], R_Z = reshape([192.05], 1, 1))

    ctrl_posicao = construirPosicaoLQR(gravidade=drone.gravidade;
        Q_X = [1232.5 0.0; 0.0 2000.0], R_X = reshape([1598.71], 1, 1),
        Q_Y = [410.2 0.0; 0.0 1270.9], R_Y = reshape([2000.00], 1, 1))
    ctrl = ControladorCascataLQR(ctrl_posicao, ctrl_altitude, ctrl_atitude)

    imprimirRelatorioLQR(ctrl)

    # ── Opção: Hélice ─────────────────────────────────────────
    tspan      = (0.0, 45.0)
    trajetoria = t -> trajetoriaHelice(t; raio=1.5, ω_giro=0.4, v_subida=0.08, t_espera=4.0)
    nome_traj  = "Hélice Ascendente"

    params = DroneTrajetoriaLQRParams(drone, ctrl, trajetoria)
    sol = resolverSistema(drone_lqr_cascata_mf!, zeros(12), tspan, params; resolucao=0.01)

    plotarTrajetoria3D(sol, trajetoria, tspan; titulo = "LQR Tracker — $nome_traj")
    plotarRespostaCompleta(sol, trajetoria; titulo = "LQR Tracker — Estados e Referências ($nome_traj)")
# =========================================================