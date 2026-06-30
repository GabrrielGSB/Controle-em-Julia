# =========================================================
# INCLUDES
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
# STRUCTS DO CONTROLADOR EM CASCATA
    struct GanhosPID
        Kp  ::Float64
        Ki  ::Float64
        Kd  ::Float64
        lim ::Float64
    end

    mutable struct ControladorCascata
        pid_Z  ::GanhosPID
        pid_X  ::GanhosPID
        pid_Y  ::GanhosPID
        pid_Φ  ::GanhosPID
        pid_θ  ::GanhosPID
        pid_Ψ  ::GanhosPID
    end

    struct DroneTrajetoriaParams
        drone      ::Drone
        ctrl       ::ControladorCascata
        trajetoria ::Function
    end

    # Dinâmica aumentada (12 estados drone + 6 integradores)
    function drone_cascata_mf!(dx, x, p::DroneTrajetoriaParams, t)
        d = p.drone;  c = p.ctrl;  g = d.gravidade

        X,  Y,  Z  = x[1], x[2], x[3]
        Vx, Vy, Vz = x[4], x[5], x[6]
        Φ,  θ,  Ψ  = x[7], x[8], x[9]
        VΦ, Vθ, VΨ = x[10], x[11], x[12]
        I_Z, I_X, I_Y = x[13], x[14], x[15]
        I_Φ, I_θ, I_Ψ = x[16], x[17], x[18]

        ref = p.trajetoria(t)
        X_r, Y_r, Z_r, Ψ_r = ref.X, ref.Y, ref.Z, ref.Ψ

        # ── Malha externa — Posição ────────────────────────────────────
        e_Z = Z_r - Z;  e_X = X_r - X;  e_Y = Y_r - Y
        Vx_i =  Vx*cos(Ψ) - Vy*sin(Ψ)
        Vy_i =  Vx*sin(Ψ) + Vy*cos(Ψ)

        U1_Δ = c.pid_Z.Kp*e_Z + c.pid_Z.Ki*I_Z + c.pid_Z.Kd*(-Vz)
        U1   = clamp(d.massa*g + U1_Δ, 0.0, 2.5*d.massa*g)

        ax_des = c.pid_X.Kp*e_X + c.pid_X.Ki*I_X + c.pid_X.Kd*(-Vx_i)
        ay_des = c.pid_Y.Kp*e_Y + c.pid_Y.Ki*I_Y + c.pid_Y.Kd*(-Vy_i)

        θ_des = clamp(( cos(Ψ)*ax_des + sin(Ψ)*ay_des)/g, -c.pid_X.lim, c.pid_X.lim)
        Φ_des = clamp(( sin(Ψ)*ax_des - cos(Ψ)*ay_des)/g, -c.pid_Y.lim, c.pid_Y.lim)

        # ── Malha interna — Atitude ────────────────────────────────────
        e_Φ = Φ_des - Φ
        e_θ = θ_des - θ
        e_Ψ = atan(sin(Ψ_r - Ψ), cos(Ψ_r - Ψ))

        U2 = clamp(c.pid_Φ.Kp*e_Φ + c.pid_Φ.Ki*I_Φ + c.pid_Φ.Kd*(-VΦ), -c.pid_Φ.lim, c.pid_Φ.lim)
        U3 = clamp(c.pid_θ.Kp*e_θ + c.pid_θ.Ki*I_θ + c.pid_θ.Kd*(-Vθ), -c.pid_θ.lim, c.pid_θ.lim)
        U4 = clamp(c.pid_Ψ.Kp*e_Ψ + c.pid_Ψ.Ki*I_Ψ + c.pid_Ψ.Kd*(-VΨ), -c.pid_Ψ.lim, c.pid_Ψ.lim)

        d.dinamica!(view(dx, 1:12), view(x, 1:12), d, t; u=[U1, U2, U3, U4])

        dx[13] = e_Z;  dx[14] = e_X;  dx[15] = e_Y
        dx[16] = e_Φ;  dx[17] = e_θ;  dx[18] = e_Ψ
    end
# =========================================================

# =========================================================
# CONFIGURAÇÃO DO ALGORITMO GENÉTICO E MÉTRICAS
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
        max_tilt   ::Float64 # Maior ângulo de Roll ou Pitch (radianos)
    end

    mutable struct Individuo
        genes    ::Vector{Float64}
        custo    ::Float64
        metricas ::Union{MetricasRastreamento, Nothing}
    end

    Individuo(genes) = Individuo(genes, Inf, nothing)
# =========================================================

# =========================================================
# LIMITES DOS GENES  (min, max) para cada um dos 18 genes
    const LIMITES_GENES = [
       0.01 30; # [1]  Kp_Z
       0   30; # [2]  Ki_Z
       0   30; # [3]  Kd_Z
       0.01 30; # [4]  Kp_X
       0   30; # [5]  Ki_X
       0   30; # [6]  Kd_X
       0.01 30; # [7]  Kp_Y
       0   30; # [8]  Ki_Y
       0   30; # [9]  Kd_Y
       0.01 30; # [10] Kp_Φ
       0   30; # [11] Ki_Φ
       0   30; # [12] Kd_Φ
       0.01 30; # [13] Kp_θ
       0   30; # [14] Ki_θ
       0   30; # [15] Kd_θ
       0.01 30; # [16] Kp_Ψ
       0   30; # [17] Ki_Ψ
       0   30; # [18] Kd_Ψ
    ]

    const LIM_ALTITUDE  = Inf
    const LIM_POSICAO   = 0.20    # ângulo desejado máximo [rad] ≈ 23°
    const LIM_TORQUE_RP = 0.60    # torque de roll/pitch [N·m]
    const LIM_TORQUE_Y  = 0.05    # torque de yaw [N·m] (Realista!)

    function decodificarGenes(g::Vector{Float64})::ControladorCascata
        return ControladorCascata(
            GanhosPID(g[1],  g[2],  g[3],  LIM_ALTITUDE),
            GanhosPID(g[4],  g[5],  g[6],  LIM_POSICAO),
            GanhosPID(g[7],  g[8],  g[9],  LIM_POSICAO),
            GanhosPID(g[10], g[11], g[12], LIM_TORQUE_RP),
            GanhosPID(g[13], g[14], g[15], LIM_TORQUE_RP),
            GanhosPID(g[16], g[17], g[18], LIM_TORQUE_Y),
        )
    end
    
    function atualizarGanhos_do_Individuo!(ctrl::ControladorCascata, ind::Individuo)
        g = ind.genes
        ctrl.pid_Z = GanhosPID(g[1],  g[2],  g[3],  LIM_ALTITUDE)
        ctrl.pid_X = GanhosPID(g[4],  g[5],  g[6],  LIM_POSICAO)
        ctrl.pid_Y = GanhosPID(g[7],  g[8],  g[9],  LIM_POSICAO)
        ctrl.pid_Φ = GanhosPID(g[10], g[11], g[12], LIM_TORQUE_RP)
        ctrl.pid_θ = GanhosPID(g[13], g[14], g[15], LIM_TORQUE_RP)
        ctrl.pid_Ψ = GanhosPID(g[16], g[17], g[18], LIM_TORQUE_Y)
    end
# =========================================================

# =========================================================
# TRAJETÓRIA EM 4 ETAPAS (WAYPOINTS SIMPLES)
    const T_SIM_AG = 20.0

    """
    Trajetória por saltos de referência (Waypoints):
    0 a 5s  : x=0.0, y=0.0, z=0.5
    5 a 10s : x=0.5, y=0.0, z=1.0
    10 a 15s: x=0.5, y=0.5, z=1.5
    15 a 20s: x=0.0, y=0.0, z=2.0
    """
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
# AVALIAÇÃO DO INDIVÍDUO E MÉTRICAS
    function avaliarIndividuo!(ind::Individuo, p::DroneTrajetoriaParams, x0, tspan)::Bool
        # 1. Injeta os genes no controlador
        atualizarGanhos_do_Individuo!(p.ctrl, ind) 

        # 2. Simula o sistema
        local sol
        try
            sol = resolverSistema(drone_cascata_mf!, x0, tspan, p; resolucao=0.05)
        catch
            ind.custo = 1e6
            return false
        end

        # 3. Punição Severa (Morte do Indivíduo) se o solver explodir ou o drone cair longe
        ultimo = sol.u[end]
        if sol.retcode != ReturnCode.Success || any(isnan, ultimo) || norm(ultimo[1:3]) > 50.0
            ind.custo = 1e6 
            return false
        end

        # 4. Avaliação Justa: CÁLCULO DO RMSE 3D AO LONGO DA ESPIRAL
        erro_quadratico_acumulado = 0.0
        max_tilt = 0.0
        
        for (i, t) in enumerate(sol.t)
            ref = p.trajetoria(t)
            X_real, Y_real, Z_real = sol.u[i][1], sol.u[i][2], sol.u[i][3]
            Roll_real, Pitch_real  = sol.u[i][7], sol.u[i][8]
            
            # Distância Euclediana ao quadrado
            erro_quadratico_acumulado += (ref.X - X_real)^2 + (ref.Y - Y_real)^2 + (ref.Z - Z_real)^2
            
            # Registra a atitude mais agressiva
            tilt_atual = max(abs(Roll_real), abs(Pitch_real))
            if tilt_atual > max_tilt
                max_tilt = tilt_atual
            end
        end
        
        rmse_3d = sqrt(erro_quadratico_acumulado / length(sol.t))

        # 5. Penalidade de Agressividade
        # Se o drone inclinar mais que ~34 graus (0.6 rad) na espiral, ele perde pontos rapidamente
        penalidade_atitude = max_tilt > 0.6 ? 50.0 * (max_tilt - 0.6) : 0.0

        # 6. Custo Final
        J = rmse_3d + penalidade_atitude
        
        ind.custo = J
        ind.metricas = MetricasRastreamento(rmse_3d, max_tilt)
        
        return true
    end
# =========================================================

# =========================================================
# OPERADORES GENÉTICOS
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
# ALGORITMO GENÉTICO PRINCIPAL
    function otimizarGanhosAG(drone::Drone; config::ConfigAG, RMSE_alvo::Float64, genes_iniciais=nothing)::Individuo
        config.semente !== nothing && Random.seed!(config.semente)
        n_pop   = config.tamanho_pop
        n_elite = max(1, round(Int, config.frac_elite * n_pop))

        println()
        println("═"^80)
        println("  AG — OTIMIZAÇÃO DE PID EM CASCATA COM ESPIRAL 3D")
        println("═"^80)
        @printf("  Meta de Parada: Erro Médio (RMSE 3D) ≤ %.3f m\n", RMSE_alvo)
        println("─"^80)

        # Inicia a população
        pop = [individuoAleatorio() for _ in 1:n_pop]
        if genes_iniciais !== nothing
            pop[1] = Individuo(copy(genes_iniciais))
        end

        # Configura o controlador e os parâmetros do solver APENAS UMA VEZ
        ctrl_mutavel = decodificarGenes(pop[1].genes)
        params_AG = DroneTrajetoriaParams(drone, ctrl_mutavel, traj_avaliacao)
        x0_AG = zeros(18)
        x0_AG[1] = 1.5 # Inicia em X=1.5 acompanhando o raio da espiral
        tspan_AG = (0.0, T_SIM_AG)

        # Avaliação Inicial
        for ind in pop
            avaliarIndividuo!(ind, params_AG, x0_AG, tspan_AG)
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
                avaliarIndividuo!(ind, params_AG, x0_AG, tspan_AG)
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

            # Critério de Parada
            if m !== nothing && m.rmse_3d ≤ RMSE_alvo && m.max_tilt ≤ 0.6
                geracao_parada = ger
                println("─"^80)
                println("  ✓ CRITÉRIO DE PARADA ATENDIDO (Drone colado na Espiral) na geração $ger")
                break
            end
        end

        println("═"^80, "\n")
        return melhor_global
    end
# =========================================================

# =========================================================
# VISUALIZAÇÃO DO AG E RESPOSTA FINAL
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

        # Plot Z
        add_trace!(fig, scatter(x=t, y=ref_Z, mode="lines", name="Ref Z", line=attr(color="black", dash="dash")), row=1, col=1)
        add_trace!(fig, scatter(x=t, y=Z, mode="lines", name="Z", line=attr(color="#1f77b4")), row=1, col=1)

        # Plot X e Y
        add_trace!(fig, scatter(x=t, y=ref_X, mode="lines", name="Ref X", line=attr(color="gray", dash="dash")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=X, mode="lines", name="X", line=attr(color="#2ca02c")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=ref_Y, mode="lines", name="Ref Y", line=attr(color="gray", dash="dot")), row=2, col=1)
        add_trace!(fig, scatter(x=t, y=Y, mode="lines", name="Y", line=attr(color="#d62728")), row=2, col=1)

        # Plot Atitude
        add_trace!(fig, scatter(x=t, y=rad2deg.(Φ), mode="lines", name="Roll Φ (°)", line=attr(color="#9467bd")), row=3, col=1)
        add_trace!(fig, scatter(x=t, y=rad2deg.(θ), mode="lines", name="Pitch θ (°)", line=attr(color="#8c564b")), row=3, col=1)

        relayout!(fig, title_text="Desempenho da Trajetória 3D Completa (Otimizada)", title_x=0.5, height=800, width=900, hovermode="x unified")
        relayout!(fig, xaxis3_title="Tempo (s)")
        display(fig)
    end
# =========================================================

# =========================================================
# EXECUÇÃO E CONFIGURAÇÃO
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

    config = ConfigAG(
        50,       # tamanho_pop
        100,      # max_geracoes
        0.20,     # frac_elite
        0.85,     # prob_crossover
        0.15,     # prob_mutacao
        0.25,     # sigma_inicial
        0.04,     # sigma_final
        3,        # torneio_k
        42        # semente
    )

    # Novo critério de Parada: Distância média (RMSE) da espiral
    alvo_rmse = 0.05 # Deseja-se que em média o drone erre no máximo 5 centímetros da linha

    genes_iniciais = [
        2.5, 0.10, 3.0,    # Z
        0.8, 0.05, 1.2,    # X
        0.8, 0.05, 1.2,    # Y
        6.0, 0.10, 1.5,    # Φ
        6.0, 0.10, 1.5,    # θ
        3.0, 0.05, 0.8     # Ψ
    ]

    # ── Execução ───────────────────────────────────────────────────────
    melhor = otimizarGanhosAG(drone; config=config, RMSE_alvo=alvo_rmse, genes_iniciais=genes_iniciais)

    # ── Validação e Gráficos ───────────────────────────────────────────
    println("A rodar validação final com a nova trajetória em 3D...\n")
    ctrl_otimo = decodificarGenes(melhor.genes)
    params_val = DroneTrajetoriaParams(drone, ctrl_otimo, traj_avaliacao)
    
    x0_val = zeros(18)

    sol_final = resolverSistema(drone_cascata_mf!, x0_val, (0.0, T_SIM_AG), params_val; resolucao=0.01)

    # Desenha os resultados no ecrã
    plotarRespostaFinal(sol_final)
# =========================================================

# =========================================================
# RELATÓRIO DOS GANHOS OTIMIZADOS
# =========================================================

g = melhor.genes

println("\n╔══════════════════════════════════════════════════════════╗")
println("║           GANHOS ÓTIMOS ENCONTRADOS PELO AG              ║")
println("╠══════════════════════════════════════════════════════════╣")
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Altitude Z", g[1],  g[2],  g[3])
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Posição X",  g[4],  g[5],  g[6])
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Posição Y",  g[7],  g[8],  g[9])
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Roll Φ",     g[10], g[11], g[12])
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Pitch θ",    g[13], g[14], g[15])
@printf("║  %-10s  Kp = %6.4f   Ki = %6.4f   Kd = %6.4f  ║\n", "Yaw Ψ",      g[16], g[17], g[18])
println("╚══════════════════════════════════════════════════════════╝\n")

println("# ── Snippet: Copie e cole no seu ficheiro principal ──────────")
@printf("ctrl_otimizado = ControladorCascata(\n")
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, Inf),\n",   g[1],  g[2],  g[3])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, %4.2f),\n",  g[4],  g[5],  g[6], LIM_POSICAO)
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, %4.2f),\n",  g[7],  g[8],  g[9], LIM_POSICAO)
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, %4.2f),\n",  g[10], g[11], g[12], LIM_TORQUE_RP)
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, %4.2f),\n",  g[13], g[14], g[15], LIM_TORQUE_RP)
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, %4.2f)\n",   g[16], g[17], g[18], LIM_TORQUE_Y)
@printf(")\n")