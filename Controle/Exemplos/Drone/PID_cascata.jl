# =========================================================
# INCLUDES
    using PlotlyJS
    using Printf
    using LinearAlgebra

    include("../../SistemasLTI/Drone.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

# =========================================================
# TRAJETÓRIAS DE REFERÊNCIA
    """
        trajetoriaHover(t; altitude, t_subida)

        Decolagem vertical suave seguida de hover estacionário.
        A subida usa perfil sigmoide para evitar descontinuidades na aceleração.
    """
    function trajetoriaHover(t; altitude=2.0, t_subida=5.0)
        τ = clamp(t / t_subida, 0.0, 1.0)
        h = altitude * τ^2 * (3 - 2τ)     # interpolação cúbica suave
        return (X=0.0, Y=0.0, Z=h, Ψ=0.0)
    end

    """
        trajetoriaHelice(t; raio, ω_giro, v_subida, t_espera)

        Decolagem vertical (t < t_espera) seguida de hélice ascendente.
    """
    function trajetoriaHelice(t; raio=1.5, ω_giro=0.4, v_subida=0.15, t_espera=4.0)
        if t < t_espera
            τ = t / t_espera
            h = τ^2 * (3 - 2τ)
            return (X=0.0, Y=0.0, Z=h, Ψ=0.0)
        end
        τ = t - t_espera
        return (X = raio * cos(ω_giro * τ),
                Y = raio * sin(ω_giro * τ),
                Z = 1.0 + v_subida * τ,
                Ψ = 0.0)
    end

    """
        trajetoriaFigura8(t; escala, T, altitude, t_subida)

        Figura-8 (lemniscata) no plano horizontal a altitude constante.
        O drone sobe suavemente antes de iniciar o percurso.
    """
    function trajetoriaFigura8(t; escala=2.0, T=20.0, altitude=2.0, t_subida=5.0)
        if t < t_subida
            τ = t / t_subida
            h = altitude * τ^2 * (3 - 2τ)
            return (X=0.0, Y=0.0, Z=h, Ψ=0.0)
        end
        τ = t - t_subida
        ω = 2π / T
        return (X = escala * sin(ω * τ),
                Y = escala * sin(ω * τ) * cos(ω * τ),
                Z = altitude,
                Ψ = 0.0)
    end

    """
        trajetoriaWaypoints(t, waypoints; T_seg)

        Sequência de waypoints com interpolação cúbica (perfil de Hermite) entre eles.

        Argumento `waypoints`: matriz (N × 4) com colunas [X, Y, Z, Ψ].
    """
    function trajetoriaWaypoints(t, waypoints; T_seg=5.0)
        n   = size(waypoints, 1)
        idx = clamp(floor(Int, t / T_seg) + 1, 1, n - 1)
        τ   = clamp((t - (idx - 1) * T_seg) / T_seg, 0.0, 1.0)
        s   = τ^2 * (3 - 2τ)                  # suavização cúbica
        p0  = waypoints[idx, :]
        p1  = waypoints[min(idx + 1, n), :]
        r   = p0 + s * (p1 - p0)
        return (X=r[1], Y=r[2], Z=r[3], Ψ=r[4])
    end
# =========================================================

# =========================================================
# STRUCTS DO CONTROLADOR EM CASCATA 
    """
        GanhosPID

        Parâmetros de um único canal PID.

        Campos:
            - Kp, Ki, Kd : Ganhos proporcional, integral e derivativo
            - lim        : Saturação simétrica da saída  |saída| ≤ lim
    """
    struct GanhosPID
        Kp  ::Float64
        Ki  ::Float64
        Kd  ::Float64
        lim ::Float64
    end

    """
        ControladorCascata

        Agrupa todos os PIDs de malha externa (posição) e interna (atitude).

        Malha externa:
            pid_Z  → Altitude       Z  → U1 (empuxo com feedforward m·g)
            pid_X  → Posição        X  → aceleração ax_des → θ_des via rotação
            pid_Y  → Posição        Y  → aceleração ay_des → Φ_des via rotação

        Malha interna:
            pid_Φ  → Roll   Φ → U2
            pid_θ  → Pitch  θ → U3
            pid_Ψ  → Yaw    Ψ → U4
    """
    struct ControladorCascata
        pid_Z  ::GanhosPID
        pid_X  ::GanhosPID
        pid_Y  ::GanhosPID
        pid_Φ  ::GanhosPID
        pid_θ  ::GanhosPID
        pid_Ψ  ::GanhosPID
    end

    """
        DroneTrajetoriaParams

        Parâmetros completos para a simulação de seguimento de trajetória.
    """
    struct DroneTrajetoriaParams
        drone      ::Drone
        ctrl       ::ControladorCascata
        trajetoria ::Function     # f(t) → NamedTuple (X, Y, Z, Ψ)
    end
# =========================================================

# =========================================================
# DINÂMICA DA MALHA FECHADA EM CASCATA
    """
        drone_cascata_mf!(dx, x, p, t)

        Equação diferencial do sistema aumentado (drone + integradores PID).

        Desacoplamento XY via rotação de Ψ:
            Para um yaw geral, a aceleração inercial desejada [ax_des, ay_des]
            é mapeada para ângulos do corpo pela inversão analítica da matriz de rotação:

                θ_des = (cos(Ψ)·ax_des + sin(Ψ)·ay_des) / g
                Φ_des = (sin(Ψ)·ax_des - cos(Ψ)·ay_des) / g

            Isso garante que os canais X e Y permaneçam desacoplados mesmo
            quando o drone está guinado.
    """
    function drone_cascata_mf!(dx, x, p::DroneTrajetoriaParams, t)
        d = p.drone
        c = p.ctrl
        g = d.gravidade

        # ── Extração dos estados do drone ─────────────────────────────
        X,  Y,  Z  = x[1], x[2], x[3]
        Vx, Vy, Vz = x[4], x[5], x[6]     # velocidades no referencial do corpo
        Φ,  θ,  Ψ  = x[7], x[8], x[9]
        VΦ, Vθ, VΨ = x[10], x[11], x[12]

        # ── Integradores ───────────────────────────────────────────────
        I_Z = x[13];  I_X = x[14];  I_Y = x[15]
        I_Φ = x[16];  I_θ = x[17];  I_Ψ = x[18]

        # ── Referência de trajetória ───────────────────────────────────
        ref             = p.trajetoria(t)
        X_r, Y_r, Z_r  = ref.X, ref.Y, ref.Z
        Ψ_r             = ref.Ψ

        # ─────────────────────────────────────────────────────────────
        # MALHA EXTERNA — CONTROLE DE POSIÇÃO
        # ─────────────────────────────────────────────────────────────

        # Erros de posição (referencial inercial)
        e_Z = Z_r - Z
        e_X = X_r - X
        e_Y = Y_r - Y

        # Velocidades no referencial inercial (rotação de Ψ)
        Vx_i =  Vx * cos(Ψ) - Vy * sin(Ψ)
        Vy_i =  Vx * sin(Ψ) + Vy * cos(Ψ)

        # PID Altitude → U1 com feedforward de gravidade
        U1_Δ = c.pid_Z.Kp * e_Z +
                c.pid_Z.Ki * I_Z +
                c.pid_Z.Kd * (-Vz)          # Vz em corpo ≈ Vz inercial (pequenos ângulos)
        U1   = clamp(d.massa * g + U1_Δ,
                    0.0, 2.5 * d.massa * g) # empuxo sempre positivo, limitado a 2.5·mg

        # PID X → aceleração inercial ax_des
        ax_des = c.pid_X.Kp * e_X +
                c.pid_X.Ki * I_X +
                c.pid_X.Kd * (-Vx_i)

        # PID Y → aceleração inercial ay_des
        ay_des = c.pid_Y.Kp * e_Y +
                c.pid_Y.Ki * I_Y +
                c.pid_Y.Kd * (-Vy_i)

        # Desacoplamento XY: rotação de Ψ para obter ângulos do corpo
        #
        #   De: X_ddot = g·(θ·cos Ψ + Φ·sin Ψ)
        #       Y_ddot = g·(θ·sin Ψ - Φ·cos Ψ)
        #
        #   Invertendo a matriz de rotação 2×2 (det = -1):
        θ_des = clamp(( cos(Ψ) * ax_des + sin(Ψ) * ay_des) / g,
                    -c.pid_X.lim, c.pid_X.lim)
        Φ_des = clamp(( sin(Ψ) * ax_des - cos(Ψ) * ay_des) / g,
                    -c.pid_Y.lim, c.pid_Y.lim)

        # ─────────────────────────────────────────────────────────────
        # MALHA INTERNA — CONTROLE DE ATITUDE
        # ─────────────────────────────────────────────────────────────

        e_Φ = Φ_des - Φ
        e_θ = θ_des - θ
        e_Ψ = atan(sin(Ψ_r - Ψ), cos(Ψ_r - Ψ))  # normalização para [-π, π]

        U2 = clamp(c.pid_Φ.Kp * e_Φ + c.pid_Φ.Ki * I_Φ + c.pid_Φ.Kd * (-VΦ),
                -c.pid_Φ.lim, c.pid_Φ.lim)

        U3 = clamp(c.pid_θ.Kp * e_θ + c.pid_θ.Ki * I_θ + c.pid_θ.Kd * (-Vθ),
                -c.pid_θ.lim, c.pid_θ.lim)

        U4 = clamp(c.pid_Ψ.Kp * e_Ψ + c.pid_Ψ.Ki * I_Ψ + c.pid_Ψ.Kd * (-VΨ),
                -c.pid_Ψ.lim, c.pid_Ψ.lim)

        # ─────────────────────────────────────────────────────────────
        # DINÂMICA DO DRONE (12 estados)
        # ─────────────────────────────────────────────────────────────
        d.dinamica!(view(dx, 1:12), view(x, 1:12), d, t; u=[U1, U2, U3, U4])

        # ─────────────────────────────────────────────────────────────
        # INTEGRADORES DOS ERROS
        # ─────────────────────────────────────────────────────────────
        dx[13] = e_Z
        dx[14] = e_X
        dx[15] = e_Y
        dx[16] = e_Φ
        dx[17] = e_θ
        dx[18] = e_Ψ
    end
# =========================================================

# =========================================================
# VISUALIZAÇÃO ESPECIALIZADA PARA TRAJETÓRIA 3D
    """
        plotarTrajetoria3D(solucao, trajetoria, tspan)

        Gera um gráfico 3D interativo comparando a trajetória do drone
        com a referência ao longo do tempo.
    """
    function plotarTrajetoria3D(solucao, trajetoria, tspan; titulo="Seguimento de Trajetória 3D")
        t_vec = solucao.t

        # Trajetória real do drone
        X_real = [u[1] for u in solucao.u]
        Y_real = [u[2] for u in solucao.u]
        Z_real = [u[3] for u in solucao.u]

        # Trajetória de referência amostrada
        refs = [trajetoria(t) for t in t_vec]
        X_ref = [r.X for r in refs]
        Y_ref = [r.Y for r in refs]
        Z_ref = [r.Z for r in refs]

        trace_real = scatter3d(
            x=X_real, y=Y_real, z=Z_real,
            mode="lines",
            name="Drone (real)",
            line=attr(color="#1f77b4", width=4)
        )
        trace_ref = scatter3d(
            x=X_ref, y=Y_ref, z=Z_ref,
            mode="lines",
            name="Referência",
            line=attr(color="#d62728", width=3, dash="dash")
        )
        trace_inicio = scatter3d(
            x=[X_real[1]], y=[Y_real[1]], z=[Z_real[1]],
            mode="markers",
            name="Início",
            marker=attr(size=6, color="green")
        )
        trace_fim = scatter3d(
            x=[X_real[end]], y=[Y_real[end]], z=[Z_real[end]],
            mode="markers",
            name="Fim",
            marker=attr(size=6, color="red", symbol="x")
        )

        layout = Layout(
            title_text=titulo, title_x=0.5,
            scene=attr(
                xaxis_title="X (m)",
                yaxis_title="Y (m)",
                zaxis_title="Z (m)"
            ),
            width=850, height=650
        )

        display(plot([trace_real, trace_ref, trace_inicio, trace_fim], layout))
    end

    """
        plotarErrosRastreamento(solucao, trajetoria)

        Plota os erros de posição e atitude ao longo do tempo.
    """
    function plotarErrosRastreamento(solucao, trajetoria; titulo="Erros de Rastreamento")
        t_vec = solucao.t

        refs   = [trajetoria(t) for t in t_vec]
        e_X    = [refs[i].X - solucao.u[i][1] for i in eachindex(t_vec)]
        e_Y    = [refs[i].Y - solucao.u[i][2] for i in eachindex(t_vec)]
        e_Z    = [refs[i].Z - solucao.u[i][3] for i in eachindex(t_vec)]
        erro_3d = [sqrt(e_X[i]^2 + e_Y[i]^2 + e_Z[i]^2) for i in eachindex(t_vec)]

        fig = make_subplots(rows=4, cols=1,
                            shared_xaxes=true,
                            vertical_spacing=0.08,
                            subplot_titles=reshape(["Erro X (m)", "Erro Y (m)",
                                                    "Erro Z (m)", "Erro Euclidiano 3D (m)"], 1, 4))

        for (linha, (dados, cor, nome)) in enumerate([
                (e_X,    "#1f77b4", "eₓ"),
                (e_Y,    "#2ca02c", "e_Y"),
                (e_Z,    "#d62728", "e_Z"),
                (erro_3d,"#9467bd", "‖e‖₂")])
            add_trace!(fig, scatter(x=t_vec, y=dados, mode="lines",
                                    name=nome, line=attr(color=cor, width=2)),
                    row=linha, col=1)
            add_trace!(fig, scatter(x=[t_vec[1], t_vec[end]], y=[0.0, 0.0],
                                    mode="lines", showlegend=false,
                                    line=attr(color="black", width=1, dash="dot")),
                    row=linha, col=1)
        end

        relayout!(fig, title_text=titulo, title_x=0.5,
                height=800, width=800, hovermode="x unified")
        relayout!(fig, xaxis4_title="Tempo (s)")
        display(fig)
    end

    """
        plotarSinaisControle(solucao, trajetoria)

        Extrai e plota U1..U4 e os ângulos desejados ao longo do tempo.
        (Exige rerodar o sistema salvando histórico — aqui recalculamos em pós-processamento.)
    """
    function plotarAtitudeEmpuxo(solucao, trajetoria)
        t_vec = solucao.t

        # Estados de posição e orientação (drone)
        X_hist  = [u[1] for u in solucao.u]
        Y_hist  = [u[2] for u in solucao.u]
        Z_hist  = [u[3] for u in solucao.u]
        Ψ_hist  = [u[9] for u in solucao.u]

        # Referências provenientes da trajetória
        refs = [trajetoria(t) for t in t_vec]
        X_ref = [r.X for r in refs]
        Y_ref = [r.Y for r in refs]
        Z_ref = [r.Z for r in refs]
        Ψ_ref = [r.Ψ for r in refs]

        fig = make_subplots(rows=4, cols=1,
                            shared_xaxes=true, vertical_spacing=0.08,
                            subplot_titles=["Posição X (m)", "Posição Y (m)",
                                            "Altitude Z (m)", "Yaw Ψ (rad)"])

        # X
        add_trace!(fig, scatter(x=t_vec, y=X_hist, name="X (drone)",
                                line=attr(color="#1f77b4", width=2)), row=1, col=1)
        add_trace!(fig, scatter(x=t_vec, y=X_ref,  name="X_ref",
                                line=attr(color="red", dash="dash", width=2)), row=1, col=1)
        # Y
        add_trace!(fig, scatter(x=t_vec, y=Y_hist, name="Y (drone)",
                                line=attr(color="#2ca02c", width=2)), row=2, col=1)
        add_trace!(fig, scatter(x=t_vec, y=Y_ref,  name="Y_ref",
                                line=attr(color="red", dash="dash", width=2)), row=2, col=1)
        # Z
        add_trace!(fig, scatter(x=t_vec, y=Z_hist, name="Z (drone)",
                                line=attr(color="#d62728", width=2)), row=3, col=1)
        add_trace!(fig, scatter(x=t_vec, y=Z_ref,  name="Z_ref",
                                line=attr(color="red", dash="dash", width=2)), row=3, col=1)
        # Ψ
        add_trace!(fig, scatter(x=t_vec, y=Ψ_hist, name="Ψ (drone)",
                                line=attr(color="#9467bd", width=2)), row=4, col=1)
        add_trace!(fig, scatter(x=t_vec, y=Ψ_ref,  name="Ψ_ref",
                                line=attr(color="red", dash="dash", width=2)), row=4, col=1)

        relayout!(fig, title_text="Estados do Drone vs Referência",
                title_x=0.5, height=900, width=800, hovermode="x unified")
        relayout!(fig, xaxis4_title="Tempo (s)")
        display(fig)
    end
# =========================================================

# =========================================================
# CONFIGURAÇÃO DO DRONE
    drone = Drone(
        massa           = 0.468,
        gravidade       = 9.81,
        Ixx             = 4.856e-3,
        Iyy             = 4.856e-3,
        Izz             = 8.801e-3,
        Ct              = 2.980e-6,
        Cl              = 1.14e-7,
        L               = 0.225,
        estadosIniciais = zeros(12)
    )
# =========================================================

# =========================================================
# GANHOS DO CONTROLADOR EM CASCATA
    ctrl = ControladorCascata(
        #           Kp     Ki    Kd     lim
        # # ── Malha externa — Posição ──────────────────────────────────
        # GanhosPID( 2.5,  0.10,  3.0,   Inf),    # Z  (altitude) → U1
        # GanhosPID( 0.8,  0.05,  1.2,   0.35),   # X  → θ_des  [rad]
        # GanhosPID( 0.8,  0.05,  1.2,   0.35),   # Y  → Φ_des  [rad]

        # # ── Malha interna — Atitude ──────────────────────────────────
        # GanhosPID( 6.0,  0.10,  1.5,   0.5),    # Φ  (roll)  → U2 [N·m]
        # GanhosPID( 6.0,  0.10,  1.5,   0.5),    # θ  (pitch) → U3 [N·m]
        # GanhosPID( 3.0,  0.05,  0.8,   0.2),    # Ψ  (yaw)   → U4 [N·m]

        # GanhosPID( 3.8587,  0.0652,  4.9918, Inf),
        # GanhosPID( 1.7021,  0.0000,  2.6208, 0.40),
        # GanhosPID( 2.2021,  0.0000,  2.4051, 0.40),
        # GanhosPID(11.4769,  0.4095,  2.4697, 0.60),
        # GanhosPID( 7.2441,  0.6217,  1.7874, 0.60),
        # GanhosPID( 2.9824,  0.2280,  1.9736, 0.25)

        GanhosPID( 7.1515,  0.0000,  3.4926, Inf),
        GanhosPID( 1.1320,  0.0000,  1.9483, 0.40),
        GanhosPID( 1.7349,  0.0000,  2.2943, 0.40),
        GanhosPID( 9.7446,  0.0000,  1.1689, 0.60),
        GanhosPID(10.0000,  0.0070,  2.6459, 0.60),
        GanhosPID( 0.0000,  0.0500,  4.0163, 0.25)
    )
# =========================================================

# =========================================================
# ESCOLHA DA TRAJETÓRIA
    # ── Opção 1: Hover estacionário (validação do controlador) ───────
    # tspan      = (0.0, 15.0)
    # trajetoria = t -> trajetoriaHover(t; altitude=2.0, t_subida=4.0)
    # nome_traj  = "Hover Estacionário (Z = 2 m)"

    # ── Opção 2: Hélice ascendente ───────────────────────────────────
    tspan      = (0.0, 40.0)
    trajetoria = t -> trajetoriaHelice(t; raio=1.5, ω_giro=0.4,
                                        v_subida=0.08, t_espera=4.0)
    nome_traj  = "Hélice Ascendente"

    # ── Opção 3: Figura-8 ────────────────────────────────────────────
    # tspan      = (0.0, 45.0)
    # trajetoria = t -> trajetoriaFigura8(t; escala=2.0, T=20.0,
    #                                        altitude=2.0, t_subida=5.0)
    # nome_traj  = "Figura-8"

    # ── Opção 4: Sequência de waypoints ──────────────────────────────
    # Waypoints  = [X  Y    Z    Ψ]
    # waypoints = [0.0  0.0  0.0  0.0;   # Solo (ponto de partida)
    #              0.0  0.0  2.0  0.0;   # Subida
    #              3.0  0.0  2.0  0.0;   # Frente
    #              3.0  3.0  3.0  0.5;   # Diagonal e subida
    #              0.0  3.0  3.0  0.5;   # Esquerda
    #              0.0  0.0  1.0  0.0]   # Retorno
    # tspan      = (0.0, 30.0)
    # trajetoria = t -> trajetoriaWaypoints(t, waypoints; T_seg=5.0)
    # nome_traj  = "Sequência de Waypoints"
# =========================================================

# =========================================================
# SIMULAÇÃO
    println("="^55)
    println("  DRONE — PID EM CASCATA")
    println("  Trajetória: $nome_traj")
    println("="^55)

    params = DroneTrajetoriaParams(drone, ctrl, trajetoria)

    # Estado inicial aumentado: [zeros(12 drone) ; zeros(6 integradores)]
    x0_aumentado = zeros(18)

    sol = resolverSistema(drone_cascata_mf!, x0_aumentado, tspan, params;
                        resolucao = 0.02)

    # ── Métricas de rastreamento ─────────────────────────────────────
    refs_final = [trajetoria(t) for t in sol.t]
    e_X_final  = [refs_final[i].X - sol.u[i][1] for i in eachindex(sol.t)]
    e_Y_final  = [refs_final[i].Y - sol.u[i][2] for i in eachindex(sol.t)]
    e_Z_final  = [refs_final[i].Z - sol.u[i][3] for i in eachindex(sol.t)]
    rmse_3d    = sqrt(mean(e_X_final.^2 + e_Y_final.^2 + e_Z_final.^2))

    @printf("\n  RMSE de rastreamento 3D : %.4f m\n",    rmse_3d)
    @printf("  Erro final em X         : %.4f m\n",    e_X_final[end])
    @printf("  Erro final em Y         : %.4f m\n",    e_Y_final[end])
    @printf("  Erro final em Z         : %.4f m\n\n",  e_Z_final[end])

    # ── Ângulos máximos atingidos ────────────────────────────────────
    Φ_max = maximum(abs.(sol.u[i][7] for i in eachindex(sol.t)))
    θ_max = maximum(abs.(sol.u[i][8] for i in eachindex(sol.t)))
    @printf("  |Φ|_max (Roll)  : %.2f° \n", rad2deg(Φ_max))
    @printf("  |θ|_max (Pitch) : %.2f° \n", rad2deg(θ_max))
    println("="^55, "\n")
# =========================================================

# =========================================================
# VISUALIZAÇÃO
    # 1. Trajetória 3D (real vs. referência)
    plotarTrajetoria3D(sol, trajetoria, tspan;
                    titulo = "Drone Cascata — $nome_traj")

    # 2. Erros de posição ao longo do tempo
    plotarErrosRastreamento(sol, trajetoria;
                            titulo = "Erros de Rastreamento — $nome_traj")

    # 3. Atitude e altitude
    plotarAtitudeEmpuxo(sol, trajetoria)
# =========================================================

# =========================================================
# ANÁLISE DE PERFORMANCE POR CANAL
    # Cria soluções fictícias com referência Z para usar analisarPerformance
    sol_z = sol   # reutiliza a mesma solução, analisando estado 3 (Z)

    m_Z = analisarPerformance(sol_z;
        referencia = trajetoria(tspan[2]).Z,  # referência final de Z
        idx_estado = 3,
        banda      = 0.05)

    imprimirRelatorio(m_Z, nome = "Canal Altitude (Z)")

    # Análise da estabilidade de Roll e Pitch (devem permanecer próximos de zero)
    m_Φ = analisarPerformance(sol;
        referencia = 0.0,
        idx_estado = 7,
        banda      = 0.10)

    m_θ = analisarPerformance(sol;
        referencia = 0.0,
        idx_estado = 8,
        banda      = 0.10)

    imprimirRelatorio(m_Φ, nome = "Estabilidade Roll Φ")
    imprimirRelatorio(m_θ, nome = "Estabilidade Pitch θ")
# =========================================================