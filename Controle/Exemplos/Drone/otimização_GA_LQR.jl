# Controle/Exemplos/Drone/otimização_GA_LQR.jl
#
# OTIMIZAÇÃO DE MATRIZES Q E R VIA ALGORITMO GENÉTICO (WAYPOINTS 3D)
# O AG busca os pesos Q/R ótimos para o LQR EM CASCATA (resolverCARE),
# exatamente como implementado em LQR_cascata.jl — SEM conversão para PID.
# Cada subsistema (Z, X, Y, Φ, θ, Ψ) é um LQR puro de 2 estados (posição/velocidade).

# =========================================================
# INCLUDES E PACOTES
    using Random
    using Statistics
    using Printf
    using PlotlyJS
    using LinearAlgebra

    include("../../SistemasLTI/Drone.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
    include("../../Ferramentas/Visualização.jl")
# =========================================================

# =========================================================
# SOLVER: EQUAÇÃO ALGÉBRICA DE RICCATI CONTÍNUA (CARE) — idêntico ao LQR_cascata.jl
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
# STRUCTS DOS SUBSISTEMAS E DO CONTROLADOR EM CASCATA — idênticos ao LQR_cascata.jl
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
# CONSTRUTORES DOS CONTROLADORES — idênticos ao LQR_cascata.jl
    function construirAtitudeLQR(Ixx::Float64, Iyy::Float64, Izz::Float64;
                                Q_roll, R_roll,
                                Q_pitch, R_pitch,
                                Q_yaw, R_yaw,
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
                                Q_Z, R_Z,
                                lim_min = 0.0, lim_max = 2.5 * massa * 9.81)
        A = [0.0  1.0; 0.0  0.0]
        B = reshape([0.0;  1.0/massa], 2, 1)
        return ControladorAltitudeLQR(SubsistemaLQR(A, B, Q_Z, R_Z; nome = "Altitude Z → U1"), massa, lim_min, lim_max)
    end

    function construirPosicaoLQR(; gravidade = 9.81,
                                Q_X, R_X,
                                Q_Y, R_Y,
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
# PARÂMETROS E DINÂMICA DA MALHA FECHADA — idênticos ao LQR_cascata.jl
#   (12 estados puros do drone, SEM estados de integrador)
    struct DroneTrajetoriaLQRParams
        drone      ::Drone
        ctrl       ::ControladorCascataLQR
        trajetoria ::Function
    end

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

        e_pos_X = [X - X_r;  Vx_i - Vx_r]
        e_pos_Y = [Y - Y_r;  Vy_i - Vy_r]

        ax_des = -(c.posicao.sub_X.K * e_pos_X)[1]
        ay_des = -(c.posicao.sub_Y.K * e_pos_Y)[1]

        θ_des = clamp(( cos(Ψ) * ax_des + sin(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)
        Φ_des = clamp(( sin(Ψ) * ax_des - cos(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)

        # ══════════════════════════════════════════════════════════════
        # NÍVEL 3 — ATITUDE
        # ══════════════════════════════════════════════════════════════
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
# CONFIGURAÇÃO DO AG E MÉTRICAS
    struct ConfigAG
        tamanho_pop    ::Int
        max_geracoes   ::Int
        frac_elite     ::Float64
        prob_crossover ::Float64
        prob_mutacao   ::Float64
        sigma_inicial  ::Float64
        sigma_final    ::Float64
        torneio_k      ::Int
        semente        ::Union{Int, Nothing}
    end

    struct MetricasRastreamento
        rmse_3d    ::Float64
        max_tilt   ::Float64
    end

    mutable struct Individuo
        genes    ::Vector{Float64}
        custo    ::Float64
        metricas ::Union{MetricasRastreamento, Nothing}
    end

    Individuo(genes) = Individuo(genes, Inf, nothing)
# =========================================================

# =========================================================
# LIMITES DOS GENES — cada canal tem 3 genes: [Q_pos, Q_vel, R]
#   (LQR puro de 2 estados, igual ao usado em construirAtitudeLQR/AltitudeLQR/PosicaoLQR)
    const LIMITES_GENES = [
       50.0  5000;  # [1]  Q_pos Z
       10.0  5000;  # [2]  Q_vel Z
        0.1  2000;  # [3]  R Z
       50.0  5000;  # [4]  Q_pos X
       10.0  5000;  # [5]  Q_vel X
        0.1  2000;  # [6]  R X
       50.0  5000;  # [7]  Q_pos Y
       10.0  5000;  # [8]  Q_vel Y
        0.1  2000;  # [9]  R Y
       50.0  5000;  # [10] Q_pos Φ
       10.0  5000;  # [11] Q_vel Φ
        0.1  2000;  # [12] R Φ
       50.0  5000;  # [13] Q_pos θ
       10.0  5000;  # [14] Q_vel θ
        0.1  2000;  # [15] R θ
       50.0  5000;  # [16] Q_pos Ψ
       10.0  5000;  # [17] Q_vel Ψ
        0.1  2000;  # [18] R Ψ
    ]

    # Constrói (de forma stateless) um ControladorCascataLQR a partir dos genes.
    # Retorna `nothing` se o Riccati falhar para alguma combinação de Q/R.
    function genesParaControladorLQR(g::Vector{Float64}, d::Drone)::Union{ControladorCascataLQR, Nothing}
        try
            ctrl_atitude = construirAtitudeLQR(d.Ixx, d.Iyy, d.Izz;
                Q_roll  = [g[10] 0.0; 0.0 g[11]], R_roll  = reshape([g[12]], 1, 1),
                Q_pitch = [g[13] 0.0; 0.0 g[14]], R_pitch = reshape([g[15]], 1, 1),
                Q_yaw   = [g[16] 0.0; 0.0 g[17]], R_yaw   = reshape([g[18]], 1, 1))

            ctrl_altitude = construirAltitudeLQR(d.massa; gravidade=d.gravidade,
                Q_Z = [g[1] 0.0; 0.0 g[2]], R_Z = reshape([g[3]], 1, 1))

            ctrl_posicao = construirPosicaoLQR(gravidade=d.gravidade;
                Q_X = [g[4] 0.0; 0.0 g[5]], R_X = reshape([g[6]], 1, 1),
                Q_Y = [g[7] 0.0; 0.0 g[8]], R_Y = reshape([g[9]], 1, 1))

            return ControladorCascataLQR(ctrl_posicao, ctrl_altitude, ctrl_atitude)
        catch
            return nothing
        end
    end
# =========================================================

# =========================================================
# TRAJETÓRIA SIMPLES (WAYPOINTS)
    const T_SIM_AG = 20.0

    function traj_avaliacao(t)
        if t < 5.0
            return (X=0.0, Y=0.0, Z=0.5, Ψ=0.0)
        elseif t < 10.0
            return (X=0.5, Y=0.0, Z=1.0, Ψ=0.0)
        elseif t < 15.0
            return (X=0.5, Y=0.5, Z=1.5, Ψ=0.0)
        else
            return (X=0.0, Y=0.0, Z=2.0, Ψ=0.0)
        end
    end
# =========================================================

# =========================================================
# AVALIAÇÃO DO INDIVÍDUO
    function avaliarIndividuo!(ind::Individuo, drone::Drone, trajetoria::Function, x0, tspan)::Bool
        # 1. Constrói o LQR em cascata a partir das matrizes Q/R do indivíduo.
        ctrl = genesParaControladorLQR(ind.genes, drone)
        if ctrl === nothing
            ind.custo = 1e6
            return false
        end

        p = DroneTrajetoriaLQRParams(drone, ctrl, trajetoria)

        # 2. Simula o sistema
        local sol
        try
            sol = resolverSistema(drone_lqr_cascata_mf!, x0, tspan, p; resolucao=0.05)
        catch
            ind.custo = 1e6
            return false
        end

        # 3. Punição severa se explodir
        ultimo = sol.u[end]
        if sol.retcode != ReturnCode.Success || any(isnan, ultimo) || norm(ultimo[1:3]) > 50.0
            ind.custo = 1e6
            return false
        end

        # 4. Cálculo do RMSE 3D (Erro Global)
        erro_quadratico_acumulado = 0.0
        max_tilt = 0.0

        for (i, t) in enumerate(sol.t)
            ref = trajetoria(t)
            X_real, Y_real, Z_real = sol.u[i][1], sol.u[i][2], sol.u[i][3]
            Roll_real, Pitch_real  = sol.u[i][7], sol.u[i][8]

            erro_quadratico_acumulado += (ref.X - X_real)^2 + (ref.Y - Y_real)^2 + (ref.Z - Z_real)^2

            tilt_atual = max(abs(Roll_real), abs(Pitch_real))
            if tilt_atual > max_tilt
                max_tilt = tilt_atual
            end
        end
        rmse_3d = sqrt(erro_quadratico_acumulado / length(sol.t))

        # 5. Penalidade de Agressividade (> 34 graus)
        penalidade_atitude = max_tilt > 0.6 ? 50.0 * (max_tilt - 0.6) : 0.0

        ind.custo = rmse_3d + penalidade_atitude
        ind.metricas = MetricasRastreamento(rmse_3d, max_tilt)
        return true
    end
# =========================================================

# =========================================================
# OPERADORES GENÉTICOS (IDÊNTICOS)
    function clipar!(genes::Vector{Float64})
        for i in eachindex(genes)
            genes[i] = clamp(genes[i], LIMITES_GENES[i, 1], LIMITES_GENES[i, 2])
        end
    end

    function individuoAleatorio()::Individuo
        genes = [LIMITES_GENES[i, 1] + rand() * (LIMITES_GENES[i, 2] - LIMITES_GENES[i, 1]) for i in 1:18]
        return Individuo(genes)
    end

    function selecaoTorneio(pop::Vector{Individuo}, k::Int)::Individuo
        concorrentes = rand(pop, k)
        return concorrentes[argmin(c.custo for c in concorrentes)]
    end

    function cruzamentoUniforme(p1::Individuo, p2::Individuo)
        g1, g2 = copy(p1.genes), copy(p2.genes)
        for i in eachindex(g1)
            if rand() < 0.5
                g1[i], g2[i] = g2[i], g1[i]
            end
        end
        return Individuo(g1), Individuo(g2)
    end

    function mutacaoGaussiana!(ind::Individuo, prob_mutacao::Float64, sigma_rel::Float64)
        for i in eachindex(ind.genes)
            if rand() < prob_mutacao
                intervalo = LIMITES_GENES[i, 2] - LIMITES_GENES[i, 1]
                ind.genes[i] += randn() * sigma_rel * intervalo
            end
        end
        clipar!(ind.genes)
    end
# =========================================================

# =========================================================
# ALGORITMO GENÉTICO LQR PRINCIPAL
    function otimizarLQR_AG(drone::Drone; config::ConfigAG, RMSE_alvo::Float64, genes_iniciais=nothing)::Individuo
        config.semente !== nothing && Random.seed!(config.semente)
        n_pop   = config.tamanho_pop
        n_elite = max(1, round(Int, config.frac_elite * n_pop))

        println()
        println("═"^80)
        println("  AG — OTIMIZAÇÃO DE MATRIZES LQR (CASCATA PURO) - WAYPOINTS 3D")
        println("═"^80)
        @printf("  Meta de Parada: Erro Médio (RMSE 3D) ≤ %.3f m\n", RMSE_alvo)
        println("─"^80)

        pop = [individuoAleatorio() for _ in 1:n_pop]
        if genes_iniciais !== nothing
            pop[1] = Individuo(copy(genes_iniciais))
        end

        x0_AG    = zeros(12)   # 12 estados puros do drone — sem integradores
        tspan_AG = (0.0, T_SIM_AG)

        for ind in pop
            avaliarIndividuo!(ind, drone, traj_avaliacao, x0_AG, tspan_AG)
        end
        sort!(pop, by=x -> x.custo)

        melhor_global = deepcopy(pop[1])
        geracao_parada = config.max_geracoes

        for ger in 1:config.max_geracoes
            σ = config.sigma_inicial + (config.sigma_final - config.sigma_inicial) * (ger - 1) / config.max_geracoes
            nova_pop = [deepcopy(pop[i]) for i in 1:n_elite]

            while length(nova_pop) < n_pop
                pai1 = selecaoTorneio(pop, config.torneio_k)
                pai2 = selecaoTorneio(pop, config.torneio_k)

                f1, f2 = rand() < config.prob_crossover ? cruzamentoUniforme(pai1, pai2) : (deepcopy(pai1), deepcopy(pai2))

                mutacaoGaussiana!(f1, config.prob_mutacao, σ)
                mutacaoGaussiana!(f2, config.prob_mutacao, σ)

                push!(nova_pop, f1)
                length(nova_pop) < n_pop && push!(nova_pop, f2)
            end

            for ind in nova_pop[n_elite+1:end]
                avaliarIndividuo!(ind, drone, traj_avaliacao, x0_AG, tspan_AG)
            end
            sort!(nova_pop, by=x -> x.custo)
            pop = nova_pop

            if pop[1].custo < melhor_global.custo
                melhor_global = deepcopy(pop[1])
            end

            m = melhor_global.metricas
            if m !== nothing
                @printf("  Ger %-3d | Custo (J): %-8.4f | RMSE 3D: %.3f m | Max Tilt: %.1f°\n",
                        ger, melhor_global.custo, m.rmse_3d, rad2deg(m.max_tilt))
            end

            if m !== nothing && m.rmse_3d ≤ RMSE_alvo && m.max_tilt ≤ 0.6
                geracao_parada = ger
                println("─"^80)
                println("  ✓ CRITÉRIO DE PARADA ATENDIDO na geração $ger")
                break
            end
        end

        println("═"^80, "\n")
        return melhor_global
    end
# =========================================================

# =========================================================
# VISUALIZAÇÃO E EXECUÇÃO
    function plotarRespostaFinal(sol)
        t    = sol.t
        X    = [u[1]  for u in sol.u]
        Y    = [u[2]  for u in sol.u]
        Z    = [u[3]  for u in sol.u]
        Φ    = [u[7]  for u in sol.u]
        θ    = [u[8]  for u in sol.u]

        ref_X = [traj_avaliacao(ti).X for ti in t]
        ref_Y = [traj_avaliacao(ti).Y for ti in t]
        ref_Z = [traj_avaliacao(ti).Z for ti in t]

        fig = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08,
                            subplot_titles=reshape(["Altitude Z (m)",
                                                    "Plano XY (m) - Navegação",
                                                    "Atitude Roll Φ e Pitch θ (rad)"], 1, 3))

        add_trace!(fig, scatter(x=t, y=ref_Z, mode="lines", name="Ref Z", line=attr(color="black", dash="dash")), row=1, col=1)
        add_trace!(fig, scatter(x=t, y=Z, mode="lines", name="Z", line=attr(color="#1f77b4")), row=1, col=1)

        add_trace!(fig, scatter(x=t, y=ref_X, mode="lines", name="Ref X", line=attr(color="gray", dash="dash")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=X, mode="lines", name="X", line=attr(color="#2ca02c")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=ref_Y, mode="lines", name="Ref Y", line=attr(color="gray", dash="dot")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=Y, mode="lines", name="Y", line=attr(color="#d62728")), row=2, col=1)

        add_trace!(fig, scatter(x=t, y=rad2deg.(Φ), mode="lines", name="Roll Φ (°)", line=attr(color="#9467bd")), row=3, col=1)
        add_trace!(fig, scatter(x=t, y=rad2deg.(θ), mode="lines", name="Pitch θ (°)", line=attr(color="#8c564b")), row=3, col=1)

        relayout!(fig, title_text="Desempenho da Trajetória (LQR Cascata Ótimo via AG)", title_x=0.5, height=800, width=900, hovermode="x unified")
        relayout!(fig, xaxis3_title="Tempo (s)")
        display(fig)
    end

    # CONFIGURAÇÃO DE EXECUÇÃO
    drone = Drone(
        massa           = 0.777,
        gravidade       = 9.81,
        Ixx             = 0.0067,
        Iyy             = 0.0059,
        Izz             = 0.0116,
        Ke              = 1.0547e-6,
        Kr              = 7.4038e-9,
        L               = 0.11,
        estadosIniciais = zeros(12)
    )

    config = ConfigAG(50, 100, 0.20, 0.85, 0.15, 0.25, 0.04, 3, 42)
    alvo_rmse = 0.05

    # Genes iniciais [Q_pos, Q_vel, R] para cada canal — semeados com os ganhos
    # já validados manualmente em LQR_cascata.jl
    genes_iniciais = [
        787.9,  796.9, 192.05,   # Z
       1232.5, 2000.0, 1598.71,  # X
        410.2, 1270.9, 2000.00,  # Y
       1028.4,  639.3,  889.30,  # Φ (roll)
       1715.3,   18.2,  290.02,  # θ (pitch)
        551.1, 1393.3,  569.36   # Ψ (yaw)
    ]

    melhor = otimizarLQR_AG(drone; config=config, RMSE_alvo=alvo_rmse, genes_iniciais=genes_iniciais)

    ctrl_otimo = genesParaControladorLQR(melhor.genes, drone)
    params_val = DroneTrajetoriaLQRParams(drone, ctrl_otimo, traj_avaliacao)
    sol_final  = resolverSistema(drone_lqr_cascata_mf!, zeros(12), (0.0, T_SIM_AG), params_val; resolucao=0.01)

    plotarRespostaFinal(sol_final)

# =========================================================

# =========================================================
# RELATÓRIO: MATRIZES Q/R E O CONTROLADOR LQR PRONTO PARA USO
# =========================================================
g = melhor.genes

println("\n╔══════════════════════════════════════════════════════════════════════════╗")
println("║                 PESOS LQR (Q e R) ENCONTRADOS PELO AG                    ║")
println("╠══════════════════════════════════════════════════════════════════════════╣")
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Altitude Z", g[1],  g[2],  g[3])
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Posição X",  g[4],  g[5],  g[6])
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Posição Y",  g[7],  g[8],  g[9])
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Roll Φ",     g[10], g[11], g[12])
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Pitch θ",    g[13], g[14], g[15])
@printf("║  %-10s  Q = diag(%7.1f, %7.1f)    R = [%7.2f]            ║\n", "Yaw Ψ",      g[16], g[17], g[18])
println("╚══════════════════════════════════════════════════════════════════════════╝\n")

println("# ── Snippet: cole diretamente em LQR_cascata.jl ──────────────────────────\n")

@printf("ctrl_atitude = construirAtitudeLQR(drone.Ixx, drone.Iyy, drone.Izz;\n")
@printf("    Q_roll  = [%.1f 0.0; 0.0 %.1f], R_roll  = reshape([%.2f], 1, 1),\n", g[10], g[11], g[12])
@printf("    Q_pitch = [%.1f 0.0; 0.0 %.1f], R_pitch = reshape([%.2f], 1, 1),\n", g[13], g[14], g[15])
@printf("    Q_yaw   = [%.1f 0.0; 0.0 %.1f], R_yaw   = reshape([%.2f], 1, 1))\n\n", g[16], g[17], g[18])

@printf("ctrl_altitude = construirAltitudeLQR(drone.massa; gravidade=drone.gravidade,\n")
@printf("    Q_Z = [%.1f 0.0; 0.0 %.1f], R_Z = reshape([%.2f], 1, 1))\n\n", g[1], g[2], g[3])

@printf("ctrl_posicao = construirPosicaoLQR(gravidade=drone.gravidade;\n")
@printf("    Q_X = [%.1f 0.0; 0.0 %.1f], R_X = reshape([%.2f], 1, 1),\n", g[4], g[5], g[6])
@printf("    Q_Y = [%.1f 0.0; 0.0 %.1f], R_Y = reshape([%.2f], 1, 1))\n", g[7], g[8], g[9])