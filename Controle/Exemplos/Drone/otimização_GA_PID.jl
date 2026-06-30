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

    struct ControladorCascata
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

    struct CriterioParada
        sobressinal_max  ::Float64
        erro_regime_max  ::Float64
    end

    struct PesosCusto
        w_os  ::Float64    # peso do sobressinal
        w_er  ::Float64    # peso do erro em regime permanente
        w_iae ::Float64    # peso do erro absoluto integral (IAE)
        w_att ::Float64    # peso da estabilidade de atitude (Roll + Pitch)
    end

    # Nova estrutura para guardar o desempenho de CADA etapa
    struct MetricasTrajetoria
        os_Z         ::Float64
        er_Z         ::Float64
        os_X         ::Float64
        er_X         ::Float64
        os_Y         ::Float64
        er_Y         ::Float64
        er_retorno_X ::Float64
        er_retorno_Y ::Float64
        iae_total    ::Float64
    end

    mutable struct Individuo
        genes    ::Vector{Float64}
        custo    ::Float64
        metricas ::Union{MetricasTrajetoria, Nothing}
    end

    Individuo(genes) = Individuo(genes, Inf, nothing)
# =========================================================

# =========================================================
# LIMITES DOS GENES  (min, max) para cada um dos 18 genes
    const LIMITES_GENES = [
       0 10; # 0.5  5.0;   # [1]  Kp_Z
       0 10; # 0.0  0.5;   # [2]  Ki_Z
       0 10; # 0.5  6.0;   # [3]  Kd_Z
       0 10; # 0.3  3.0;   # [4]  Kp_X
       0 10; # 0.0  0.3;   # [5]  Ki_X
       0 10; # 0.3  4.0;   # [6]  Kd_X
       0 10; # 0.3  3.0;   # [7]  Kp_Y
       0 10; # 0.0  0.3;   # [8]  Ki_Y
       0 10; # 0.3  4.0;   # [9]  Kd_Y
       0 10; # 2.0  15.0;  # [10] Kp_Φ
       0 10; # 0.0  1.0;   # [11] Ki_Φ
       0 10; # 0.5  5.0;   # [12] Kd_Φ
       0 10; # 2.0  15.0;  # [13] Kp_θ
       0 10; # 0.0  1.0;   # [14] Ki_θ
       0 10; # 0.5  5.0;   # [15] Kd_θ
       0 10; # 1.0  8.0;   # [16] Kp_Ψ
       0 10; # 0.0  0.5;   # [17] Ki_Ψ
       0 10;# 0.3  3.0;   # [18] Kd_Ψ
    ]

    const LIM_ALTITUDE  = Inf
    const LIM_POSICAO   = 0.40    # ângulo desejado máximo [rad] ≈ 23°
    const LIM_TORQUE_RP = 0.60    # torque de roll/pitch [N·m]
    const LIM_TORQUE_Y  = 0.25    # torque de yaw [N·m]

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
# =========================================================

# =========================================================
# TRAJETÓRIA EM 4 ETAPAS 
    const T_SIM_AG = 20.0

    """
    Trajetória faseada:
    0 a 5s  : Sobe para Z = 1m
    5 a 10s : Anda para X = 1m
    10 a 15s: Anda para Y = 1m
    15 a 20s: Retorna diagonalmente para X = 0, Y = 0
    """
    function traj_avaliacao(t)
        if t < 5.0
            return (X=0.0, Y=0.0, Z=1.0, Ψ=0.0)
        elseif t < 10.0
            return (X=1.0, Y=0.0, Z=1.0, Ψ=0.0)
        elseif t < 15.0
            return (X=1.0, Y=1.0, Z=1.0, Ψ=0.0)
        else
            return (X=0.0, Y=0.0, Z=1.0, Ψ=0.0)
        end
    end
# =========================================================

# =========================================================
# AVALIAÇÃO DO INDIVÍDUO E MÉTRICAS POR FASE
    """
        Extrai o sobressinal (%) e o erro de regime (m) de uma variável específica, 
        apenas dentro do intervalo de tempo [t_inicio, t_fim].
    """
    function extrair_metricas_fase(sol, t_inicio, t_fim, idx_estado, ref_local)
        indices = findall(t -> t_inicio <= t <= t_fim, sol.t)
        if isempty(indices) return 0.0, Inf end
        
        valores = [sol.u[i][idx_estado] for i in indices]
        val_inicial = valores[1]
        val_final = valores[end]
        degrau = ref_local - val_inicial
        
        if abs(degrau) > 1e-3
            if degrau > 0
                pico = maximum(valores)
                os = max(0.0, (pico - ref_local) / degrau) * 100
            else
                pico = minimum(valores)
                os = max(0.0, (ref_local - pico) / abs(degrau)) * 100
            end
        else
            os = 0.0
        end
        
        er = abs(val_final - ref_local)
        return os, er
    end

    function avaliarIndividuo!(ind::Individuo, drone::Drone, pesos::PesosCusto)::Bool
        ctrl   = decodificarGenes(ind.genes)
        params = DroneTrajetoriaParams(drone, ctrl, traj_avaliacao)

        local sol
        try
            sol = resolverSistema(drone_cascata_mf!, zeros(18),
                                (0.0, T_SIM_AG), params; resolucao=0.05)
        catch
            ind.custo = Inf
            return false
        end

        ultimo = sol.u[end]
        if any(isnan, ultimo) || any(isinf, ultimo) || norm(ultimo[1:3]) > 50.0
            ind.custo = Inf
            return false
        end

        Φ_max = maximum(abs(sol.u[i][7]) for i in eachindex(sol.t))
        θ_max = maximum(abs(sol.u[i][8]) for i in eachindex(sol.t))
        if Φ_max > deg2rad(60.0) || θ_max > deg2rad(60.0)
            ind.custo = Inf
            return false
        end

        # ── Avaliação de cada etapa da Trajetória ────────────────────────
        # 1. Subida Z (0 a 5s)
        os_Z, er_Z = extrair_metricas_fase(sol, 0.0, 5.0, 3, 1.0)
        
        # 2. Movimento X (5 a 10s)
        os_X, er_X = extrair_metricas_fase(sol, 5.0, 10.0, 1, 1.0)
        
        # 3. Movimento Y (10 a 15s)
        os_Y, er_Y = extrair_metricas_fase(sol, 10.0, 15.0, 2, 1.0)
        
        # 4. Retorno Diagonal XY (15 a 20s)
        _, er_ret_X = extrair_metricas_fase(sol, 15.0, 20.0, 1, 0.0)
        _, er_ret_Y = extrair_metricas_fase(sol, 15.0, 20.0, 2, 0.0)

        # IAE Total (Área do erro para rastreamento suave)
        iae_total = sum(
            (abs(traj_avaliacao(sol.t[i]).X - sol.u[i][1]) + 
            abs(traj_avaliacao(sol.t[i]).Y - sol.u[i][2]) + 
            abs(traj_avaliacao(sol.t[i]).Z - sol.u[i][3])) * 0.05 
            for i in eachindex(sol.t)
        )

        # RMS de Atitude (Penaliza oscilações de Roll e Pitch)
        Φ_rms = sqrt(mean(sol.u[i][7]^2 for i in eachindex(sol.t)))
        θ_rms = sqrt(mean(sol.u[i][8]^2 for i in eachindex(sol.t)))
        att_rms_deg = rad2deg(Φ_rms + θ_rms)

        # ── Função de Custo Agregada ──────────────────────────────────────
        # A penalidade abrange a performance geral de todas as fases.
        custo_os = os_Z^2 + os_X^2 + os_Y^2
        custo_er = er_Z^2 + er_X^2 + er_Y^2 + er_ret_X^2 + er_ret_Y^2

        ind.custo = (pesos.w_os * custo_os + 
                    pesos.w_er * custo_er + 
                    pesos.w_iae * iae_total + 
                    pesos.w_att * att_rms_deg)

        ind.metricas = MetricasTrajetoria(os_Z, er_Z, os_X, er_X, os_Y, er_Y, er_ret_X, er_ret_Y, iae_total)
        
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
# VERIFICAÇÃO DO CRITÉRIO DE PARADA
    """
        A simulação só para precocemente se a performance for decente
        em TODAS as fases da trajetória.
    """
    function criterioAtendido(ind::Individuo, crit::CriterioParada)::Bool
        ind.metricas === nothing && return false
        m = ind.metricas

        os_ok = (m.os_Z ≤ crit.sobressinal_max) && (m.os_X ≤ crit.sobressinal_max) && (m.os_Y ≤ crit.sobressinal_max)
        er_ok = (m.er_Z ≤ crit.erro_regime_max) && (m.er_X ≤ crit.erro_regime_max) && (m.er_Y ≤ crit.erro_regime_max) && 
                (m.er_retorno_X ≤ crit.erro_regime_max) && (m.er_retorno_Y ≤ crit.erro_regime_max)

        return os_ok && er_ok
    end
# =========================================================

# =========================================================
# VISUALIZAÇÃO DO AG E RESPOSTA FINAL
    function plotarRespostaFinal(sol, criterio::CriterioParada)
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
# ALGORITMO GENÉTICO PRINCIPAL
    function otimizarGanhosAG(drone::Drone; config::ConfigAG, criterio::CriterioParada, pesos::PesosCusto, genes_iniciais=nothing)::Individuo
        config.semente !== nothing && Random.seed!(config.semente)
        n_pop   = config.tamanho_pop
        n_elite = max(1, round(Int, config.frac_elite * n_pop))

        println()
        println("═"^80)
        println("  AG — OTIMIZAÇÃO DE PID EM CASCATA COM MANOBRA 3D")
        println("═"^80)
        @printf("  Metas para TODAS as Fases: OS ≤ %.1f%%  |  Erro Regime ≤ %.3f m\n", 
                criterio.sobressinal_max, criterio.erro_regime_max)
        println("─"^80)

        pop = [individuoAleatorio() for _ in 1:n_pop]
        if genes_iniciais !== nothing
            pop[1] = Individuo(copy(genes_iniciais))
        end

        for ind in pop
            avaliarIndividuo!(ind, drone, pesos)
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
                avaliarIndividuo!(ind, drone, pesos)
            end
            sort!(nova_pop, by=x -> x.custo)
            pop = nova_pop

            if pop[1].custo < melhor_global.custo
                melhor_global = deepcopy(pop[1])
            end

            m = melhor_global.metricas
            @printf("  Ger %-3d | Custo: %-8.2f | OS(Z,X,Y): %4.1f%%, %4.1f%%, %4.1f%% | IAE: %.2f\n",
                    ger, melhor_global.custo, m.os_Z, m.os_X, m.os_Y, m.iae_total)

            if criterioAtendido(melhor_global, criterio)
                geracao_parada = ger
                println("─"^80)
                println("  ✓ CRITÉRIO DE PARADA ATENDIDO (TODAS AS FASES DECENTES) na geração $ger")
                break
            end
        end

        println("═"^80, "\n")
        return melhor_global
    end
# =========================================================

# =========================================================
# EXECUÇÃO E CONFIGURAÇÃO
    drone = Drone(
        massa = 0.468, gravidade = 9.81, Ixx = 4.856e-3, Iyy = 4.856e-3, Izz = 8.801e-3,
        Ct = 2.980e-6, Cl = 1.14e-7, L = 0.225, estadosIniciais = zeros(12)
    )

    config = ConfigAG(
        50,       # tamanho_pop
        100,      # max_geracoes
        0.20,     # frac_elite
        0.85,     # prob_crossover
        0.15,     # prob_mutacao
        0.25,     # sigma_inicial (ligeiramente maior para explorar o espaço 3D)
        0.04,     # sigma_final
        3,        # torneio_k
        42        # semente
    )

    # Critério rígido: queremos overshoot menor que 10% e erro menor que 5cm em TODAS as transições
    criterio = CriterioParada(10.0, 0.05)

    pesos = PesosCusto(
        2.0,    # w_os  : penaliza fortemente sobressinal
        10.0,   # w_er  : penaliza fortemente erro em regime (em todos os eixos)
        0.5,    # w_iae : área do erro
        2.5     # w_att : penaliza balanço extremo de Roll e Pitch nas curvas
    )

    genes_iniciais = [
        2.5, 0.10, 3.0,    # Z
        0.8, 0.05, 1.2,    # X
        0.8, 0.05, 1.2,    # Y
        6.0, 0.10, 1.5,    # Φ
        6.0, 0.10, 1.5,    # θ
        3.0, 0.05, 0.8     # Ψ
    ]

    # ── Execução ───────────────────────────────────────────────────────
    melhor = otimizarGanhosAG(drone; config=config, criterio=criterio, pesos=pesos, genes_iniciais=genes_iniciais)

    # ── Validação e Gráficos ───────────────────────────────────────────
    println("A rodar validação final com a nova trajetória em 3D...\n")
    ctrl_otimo = decodificarGenes(melhor.genes)
    params_val = DroneTrajetoriaParams(drone, ctrl_otimo, traj_avaliacao)

    sol_final = resolverSistema(drone_cascata_mf!, zeros(18), (0.0, T_SIM_AG), params_val; resolucao=0.01)

    # Desenha os resultados no ecrã
    plotarRespostaFinal(sol_final, criterio)
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

# Snippet opcional para copiar e colar diretamente noutros ficheiros do seu projeto:
println("# ── Snippet: Copie e cole no seu ficheiro principal ──────────")
@printf("ctrl_otimizado = ControladorCascata(\n")
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, Inf),\n",   g[1],  g[2],  g[3])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, 0.40),\n",  g[4],  g[5],  g[6])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, 0.40),\n",  g[7],  g[8],  g[9])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, 0.60),\n",  g[10], g[11], g[12])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, 0.60),\n",  g[13], g[14], g[15])
@printf("    GanhosPID(%7.4f, %7.4f, %7.4f, 0.25)\n",   g[16], g[17], g[18])
@printf(")\n")