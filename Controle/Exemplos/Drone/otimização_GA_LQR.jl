# Controle/Exemplos/Drone/AG_Otimizacao_LQR_3D.jl
#
# OTIMIZAÇÃO DE MATRIZES Q E R (LQR EM CASCATA) VIA ALGORITMO GENÉTICO (MANOBRA 3D)
# COM FORMULAÇÃO LQR TRACKER (RASTREADOR DE REFERÊNCIA)

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
# SOLVER CARE E STRUCTS DO LQR
# =========================================================
function resolverCARE(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64})
    n = size(A, 1)
    Rinv = inv(R)
    H = [ A           -B*Rinv*B';
         -Q           -A'       ]
    vals, vecs = eigen(H)
    ordem = sortperm(real.(vals))
    estavel = ordem[1:n]

    if any(real.(vals[estavel]) .≥ -1e-9)
        error("CARE sem autovalores estritamente estáveis.")
    end

    V = vecs[:, estavel]
    V1, V2 = V[1:n, :], V[n+1:end, :]
    P = real(V2 / V1)
    P = Matrix(Symmetric((P + P') / 2))
    K = Rinv * B' * P
    polos = eigvals(A - B * K)
    return K, P, polos
end

struct SubsistemaLQR
    A::Matrix{Float64}; B::Matrix{Float64}
    Q::Matrix{Float64}; R::Matrix{Float64}
    K::Matrix{Float64}; P::Matrix{Float64}
    polos::Vector{ComplexF64}; nome::String
end

function SubsistemaLQR(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}; nome::String = "Sub")
    K, P, polos = resolverCARE(A, B, Q, R)
    return SubsistemaLQR(A, B, Q, R, K, P, polos, nome)
end

struct ControladorAtitudeLQR
    roll::SubsistemaLQR; pitch::SubsistemaLQR; yaw::SubsistemaLQR
    lim_U2::Float64; lim_U3::Float64; lim_U4::Float64
end

struct ControladorAltitudeLQR
    sub::SubsistemaLQR
    massa::Float64; lim_min::Float64; lim_max::Float64
end

struct ControladorPosicaoLQR
    sub_X::SubsistemaLQR; sub_Y::SubsistemaLQR
    lim_angulo::Float64; gravidade::Float64
end

struct ControladorCascataLQR
    posicao::ControladorPosicaoLQR
    altitude::ControladorAltitudeLQR
    atitude::ControladorAtitudeLQR
end

struct DroneTrajetoriaLQRParams
    drone::Drone
    ctrl::ControladorCascataLQR
    trajetoria::Function
end

# Dinâmica do LQR TRACKER (Atua no Erro de Estado e de Velocidade)
function drone_lqr_cascata_mf!(dx, x, p::DroneTrajetoriaLQRParams, t)
    d = p.drone; c = p.ctrl; g = d.gravidade

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

    # Nível 2 - Altitude
    e_alt = [Z - Z_r; Vz - Vz_r]
    δU1   = -(c.altitude.sub.K * e_alt)[1]
    U1    = clamp(c.altitude.massa * g + δU1, c.altitude.lim_min, c.altitude.lim_max)

    # Nível 1 - Posição
    Vx_i =  Vx * cos(Ψ) - Vy * sin(Ψ)
    Vy_i =  Vx * sin(Ψ) + Vy * cos(Ψ)
    e_pos_X = [X - X_r; Vx_i - Vx_r]
    e_pos_Y = [Y - Y_r; Vy_i - Vy_r]

    ax_des = -(c.posicao.sub_X.K * e_pos_X)[1]
    ay_des = -(c.posicao.sub_Y.K * e_pos_Y)[1]

    θ_des = clamp(( cos(Ψ) * ax_des + sin(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)
    Φ_des = clamp(( sin(Ψ) * ax_des - cos(Ψ) * ay_des) / g, -c.posicao.lim_angulo, c.posicao.lim_angulo)

    # Nível 3 - Atitude
    e_roll  = [Φ - Φ_des; VΦ - 0.0]
    e_pitch = [θ - θ_des; Vθ - 0.0]
    e_yaw   = [atan(sin(Ψ - Ψ_r), cos(Ψ - Ψ_r)); VΨ - VΨ_r]

    U2 = clamp(-(c.atitude.roll.K  * e_roll )[1], -c.atitude.lim_U2, c.atitude.lim_U2)
    U3 = clamp(-(c.atitude.pitch.K * e_pitch)[1], -c.atitude.lim_U3, c.atitude.lim_U3)
    U4 = clamp(-(c.atitude.yaw.K   * e_yaw  )[1], -c.atitude.lim_U4, c.atitude.lim_U4)

    d.dinamica!(dx, x, d, t; u = [U1, U2, U3, U4])
end

# =========================================================
# CONFIGURAÇÃO DO ALGORITMO GENÉTICO
# =========================================================

struct ConfigAG
    tamanho_pop::Int; max_geracoes::Int; frac_elite::Float64
    prob_crossover::Float64; prob_mutacao::Float64
    sigma_inicial::Float64; sigma_final::Float64
    torneio_k::Int; semente::Union{Int, Nothing}
end

struct CriterioParada
    sobressinal_max::Float64
    erro_regime_max::Float64
end

struct PesosCusto
    w_os::Float64; w_er::Float64; w_iae::Float64; w_att::Float64
end

struct MetricasTrajetoria
    os_Z::Float64; er_Z::Float64; os_X::Float64; er_X::Float64
    os_Y::Float64; er_Y::Float64; er_retorno_X::Float64; er_retorno_Y::Float64
    iae_total::Float64
end

mutable struct Individuo
    genes::Vector{Float64}
    custo::Float64
    metricas::Union{MetricasTrajetoria, Nothing}
end
Individuo(genes) = Individuo(genes, Inf, nothing)

# =========================================================
# LIMITES DOS GENES (Q11, Q22, R para as 6 malhas)
# =========================================================
if !@isdefined(LIMITES_GENES)
    const LIMITES_GENES = [
        # Z (Altitude)
        0.0 2000.0; # [1]  Q_Z_pos
        0.0 2000.0; # [2]  Q_Z_vel
        0.0 2000.0; # [3]  R_Z
        # X (Posição)
        0.0 2000.0; # [4]  Q_X_pos
        0.0 2000.0; # [5]  Q_X_vel
        0.0 2000.0; # [6]  R_X
        # Y (Posição)
        0.0 2000.0; # [7]  Q_Y_pos
        0.0 2000.0; # [8]  Q_Y_vel
        0.0 2000.0; # [9]  R_Y
        # Roll (Φ)
        0.0 2000.0; # [10] Q_Φ_ang
        0.0 2000.0; # [11] Q_Φ_taxa
        0.0 2000.0; # [12] R_Φ
        # Pitch (θ)
        0.0 2000.0; # [13] Q_θ_ang
        0.0 2000.0; # [14] Q_θ_taxa
        0.0 2000.0; # [15] R_θ
        # Yaw (Ψ)
        0.0 2000.0; # [16] Q_Ψ_ang
        0.0 2000.0; # [17] Q_Ψ_taxa
        0.0 2000.0; # [18] R_Ψ
    ]
end

function decodificarGenesLQR(g::Vector{Float64}, d::Drone)::Union{ControladorCascataLQR, Nothing}
    A_di = [0.0 1.0; 0.0 0.0]
    
    try
        # Altitude
        B_Z = reshape([0.0; 1.0/d.massa], 2, 1)
        sub_Z = SubsistemaLQR(A_di, B_Z, [g[1] 0.0; 0.0 g[2]], reshape([g[3]], 1, 1))
        ctrl_Z = ControladorAltitudeLQR(sub_Z, d.massa, 0.0, 2.5 * d.massa * d.gravidade)

        # Posicao
        B_XY = reshape([0.0; 1.0], 2, 1)
        sub_X = SubsistemaLQR(A_di, B_XY, [g[4] 0.0; 0.0 g[5]], reshape([g[6]], 1, 1))
        sub_Y = SubsistemaLQR(A_di, B_XY, [g[7] 0.0; 0.0 g[8]], reshape([g[9]], 1, 1))
        ctrl_pos = ControladorPosicaoLQR(sub_X, sub_Y, 0.35, d.gravidade)

        # Atitude
        B_roll = reshape([0.0; 1.0/d.Ixx], 2, 1)
        B_pitch = reshape([0.0; 1.0/d.Iyy], 2, 1)
        B_yaw = reshape([0.0; 1.0/d.Izz], 2, 1)
        
        roll  = SubsistemaLQR(A_di, B_roll, [g[10] 0.0; 0.0 g[11]], reshape([g[12]], 1, 1))
        pitch = SubsistemaLQR(A_di, B_pitch, [g[13] 0.0; 0.0 g[14]], reshape([g[15]], 1, 1))
        yaw   = SubsistemaLQR(A_di, B_yaw, [g[16] 0.0; 0.0 g[17]], reshape([g[18]], 1, 1))
        ctrl_att = ControladorAtitudeLQR(roll, pitch, yaw, 0.60, 0.60, 0.25)

        return ControladorCascataLQR(ctrl_pos, ctrl_Z, ctrl_att)
    catch
        return nothing # Falha no Riccati (CARE)
    end
end

# =========================================================
# TRAJETÓRIA E AVALIAÇÃO
# =========================================================
const T_SIM_AG = 20.0

function traj_avaliacao(t)
    # Como são degraus puros, a velocidade de referência é 0 em todos os eixos.
    if t < 5.0;      return (X=0.0, Y=0.0, Z=1.0, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=0.0, VΨ=0.0)
    elseif t < 10.0; return (X=1.0, Y=0.0, Z=1.0, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=0.0, VΨ=0.0)
    elseif t < 15.0; return (X=1.0, Y=1.0, Z=1.0, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=0.0, VΨ=0.0)
    else;            return (X=0.0, Y=0.0, Z=1.0, Ψ=0.0, Vx=0.0, Vy=0.0, Vz=0.0, VΨ=0.0)
    end
end

function extrair_metricas_fase(sol, t_inicio, t_fim, idx_estado, ref_local)
    indices = findall(t -> t_inicio <= t <= t_fim, sol.t)
    if isempty(indices) return 0.0, Inf end
    valores = [sol.u[i][idx_estado] for i in indices]
    degrau = ref_local - valores[1]
    
    if abs(degrau) > 1e-3
        os = degrau > 0 ? max(0.0, (maximum(valores) - ref_local) / degrau) * 100 : 
                          max(0.0, (ref_local - minimum(valores)) / abs(degrau)) * 100
    else
        os = 0.0
    end
    return os, abs(valores[end] - ref_local)
end

function avaliarIndividuo!(ind::Individuo, drone::Drone, pesos::PesosCusto)::Bool
    ctrl = decodificarGenesLQR(ind.genes, drone)
    
    # 1. Trava do Solver (Falha do Riccati)
    if ctrl === nothing
        ind.custo = Inf
        return false
    end

    # 2. Trava da Arquitetura em Cascata (Separação de Escala de Tempo)
    ωn_att = abs(ctrl.atitude.roll.polos[1])
    ωn_pos = abs(ctrl.posicao.sub_X.polos[1])
    if ωn_pos * 4 > ωn_att  # A atitude tem que ser PELO MENOS 4x mais rápida que a posição
        ind.custo = Inf
        return false
    end

    params = DroneTrajetoriaLQRParams(drone, ctrl, traj_avaliacao)
    local sol
    try
        sol = resolverSistema(drone_lqr_cascata_mf!, zeros(12), (0.0, T_SIM_AG), params; resolucao=0.05)
    catch
        ind.custo = Inf; return false
    end

    ultimo = sol.u[end]
    if any(isnan, ultimo) || any(isinf, ultimo) || norm(ultimo[1:3]) > 50.0
        ind.custo = Inf; return false
    end

    Φ_max = maximum(abs(sol.u[i][7]) for i in eachindex(sol.t))
    θ_max = maximum(abs(sol.u[i][8]) for i in eachindex(sol.t))
    if Φ_max > deg2rad(60.0) || θ_max > deg2rad(60.0)
        ind.custo = Inf; return false
    end

    # Métricas
    os_Z, er_Z = extrair_metricas_fase(sol, 0.0, 5.0, 3, 1.0)
    os_X, er_X = extrair_metricas_fase(sol, 5.0, 10.0, 1, 1.0)
    os_Y, er_Y = extrair_metricas_fase(sol, 10.0, 15.0, 2, 1.0)
    _, er_ret_X = extrair_metricas_fase(sol, 15.0, 20.0, 1, 0.0)
    _, er_ret_Y = extrair_metricas_fase(sol, 15.0, 20.0, 2, 0.0)

    iae_total = sum((abs(traj_avaliacao(sol.t[i]).X - sol.u[i][1]) + 
                     abs(traj_avaliacao(sol.t[i]).Y - sol.u[i][2]) + 
                     abs(traj_avaliacao(sol.t[i]).Z - sol.u[i][3])) * 0.05 for i in eachindex(sol.t))

    Φ_rms = sqrt(mean(sol.u[i][7]^2 for i in eachindex(sol.t)))
    θ_rms = sqrt(mean(sol.u[i][8]^2 for i in eachindex(sol.t)))
    att_rms_deg = rad2deg(Φ_rms + θ_rms)

    custo_os = os_Z^2 + os_X^2 + os_Y^2
    custo_er = er_Z^2 + er_X^2 + er_Y^2 + er_ret_X^2 + er_ret_Y^2

    ind.custo = (pesos.w_os * custo_os + pesos.w_er * custo_er + pesos.w_iae * iae_total + pesos.w_att * att_rms_deg)
    ind.metricas = MetricasTrajetoria(os_Z, er_Z, os_X, er_X, os_Y, er_Y, er_ret_X, er_ret_Y, iae_total)
    
    return true
end

# =========================================================
# OPERADORES GENÉTICOS E MAIN
# =========================================================
function clipar!(genes::Vector{Float64})
    for i in eachindex(genes)
        genes[i] = clamp(genes[i], LIMITES_GENES[i, 1], LIMITES_GENES[i, 2])
    end
end

function individuoAleatorio()::Individuo
    return Individuo([LIMITES_GENES[i, 1] + rand() * (LIMITES_GENES[i, 2] - LIMITES_GENES[i, 1]) for i in 1:18])
end

function selecaoTorneio(pop::Vector{Individuo}, k::Int)::Individuo
    return rand(pop, k) |> (c -> c[argmin([x.custo for x in c])])
end

function cruzamentoUniforme(p1::Individuo, p2::Individuo)
    g1, g2 = copy(p1.genes), copy(p2.genes)
    for i in eachindex(g1)
        if rand() < 0.5; g1[i], g2[i] = g2[i], g1[i]; end
    end
    return Individuo(g1), Individuo(g2)
end

function mutacaoGaussiana!(ind::Individuo, prob_mutacao::Float64, sigma_rel::Float64)
    for i in eachindex(ind.genes)
        if rand() < prob_mutacao
            ind.genes[i] += randn() * sigma_rel * (LIMITES_GENES[i, 2] - LIMITES_GENES[i, 1])
        end
    end
    clipar!(ind.genes)
end

function criterioAtendido(ind::Individuo, crit::CriterioParada)::Bool
    ind.metricas === nothing && return false
    m = ind.metricas
    os_ok = (m.os_Z ≤ crit.sobressinal_max) && (m.os_X ≤ crit.sobressinal_max) && (m.os_Y ≤ crit.sobressinal_max)
    er_ok = (m.er_Z ≤ crit.erro_regime_max) && (m.er_X ≤ crit.erro_regime_max) && (m.er_Y ≤ crit.erro_regime_max) && 
            (m.er_retorno_X ≤ crit.erro_regime_max) && (m.er_retorno_Y ≤ crit.erro_regime_max)
    return os_ok && er_ok
end

function plotarRespostaFinal(sol, criterio::CriterioParada)
    t = sol.t
    X = [u[1] for u in sol.u]; Y = [u[2] for u in sol.u]; Z = [u[3] for u in sol.u]
    Φ = [u[7] for u in sol.u]; θ = [u[8] for u in sol.u]
    
    fig = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08,
                        subplot_titles=reshape(["Altitude Z", "Plano XY - Navegação", "Atitude Roll e Pitch"], 1, 3))

    add_trace!(fig, scatter(x=t, y=[traj_avaliacao(ti).Z for ti in t], name="Ref Z", line=attr(color="black", dash="dash")), row=1, col=1)
    add_trace!(fig, scatter(x=t, y=Z, name="Z", line=attr(color="#1f77b4")), row=1, col=1)

    add_trace!(fig, scatter(x=t, y=[traj_avaliacao(ti).X for ti in t], name="Ref X", line=attr(color="gray", dash="dash")), row=2, col=1)
    add_trace!(fig, scatter(x=t, y=X, name="X", line=attr(color="#2ca02c")), row=2, col=1)
    add_trace!(fig, scatter(x=t, y=[traj_avaliacao(ti).Y for ti in t], name="Ref Y", line=attr(color="gray", dash="dot")), row=2, col=1)
    add_trace!(fig, scatter(x=t, y=Y, name="Y", line=attr(color="#d62728")), row=2, col=1)

    add_trace!(fig, scatter(x=t, y=rad2deg.(Φ), name="Roll Φ (°)", line=attr(color="#9467bd")), row=3, col=1)
    add_trace!(fig, scatter(x=t, y=rad2deg.(θ), name="Pitch θ (°)", line=attr(color="#8c564b")), row=3, col=1)

    relayout!(fig, title_text="LQR Otimizado - Trajetória 3D", title_x=0.5, height=800, width=900, hovermode="x unified")
    display(fig)
end

function otimizarLQRAg(drone::Drone; config::ConfigAG, criterio::CriterioParada, pesos::PesosCusto, genes_iniciais=nothing)::Individuo
    config.semente !== nothing && Random.seed!(config.semente)
    n_pop = config.tamanho_pop
    n_elite = max(1, round(Int, config.frac_elite * n_pop))

    println("\n" * "═"^80)
    println("  AG — OTIMIZAÇÃO DE MATRIZES Q E R (LQR CASCATA TRACKER)")
    println("═"^80)

    pop = [individuoAleatorio() for _ in 1:n_pop]
    if genes_iniciais !== nothing; pop[1] = Individuo(copy(genes_iniciais)); end

    for ind in pop; avaliarIndividuo!(ind, drone, pesos); end
    sort!(pop, by=x -> x.custo)

    melhor_global = deepcopy(pop[1])

    for ger in 1:config.max_geracoes
        σ = config.sigma_inicial + (config.sigma_final - config.sigma_inicial) * (ger - 1) / config.max_geracoes
        nova_pop = [deepcopy(pop[i]) for i in 1:n_elite]

        while length(nova_pop) < n_pop
            pai1 = selecaoTorneio(pop, config.torneio_k)
            pai2 = selecaoTorneio(pop, config.torneio_k)
            f1, f2 = rand() < config.prob_crossover ? cruzamentoUniforme(pai1, pai2) : (deepcopy(pai1), deepcopy(pai2))
            mutacaoGaussiana!(f1, config.prob_mutacao, σ); mutacaoGaussiana!(f2, config.prob_mutacao, σ)
            push!(nova_pop, f1); length(nova_pop) < n_pop && push!(nova_pop, f2)
        end

        for ind in nova_pop[n_elite+1:end]; avaliarIndividuo!(ind, drone, pesos); end
        sort!(nova_pop, by=x -> x.custo)
        pop = nova_pop

        if pop[1].custo < melhor_global.custo
            melhor_global = deepcopy(pop[1])
        end

        m = melhor_global.metricas
        @printf("  Ger %-3d | Custo: %-8.2f | OS(Z,X,Y): %4.1f%%, %4.1f%%, %4.1f%% | IAE: %.2f\n",
                ger, melhor_global.custo, m.os_Z, m.os_X, m.os_Y, m.iae_total)

        if criterioAtendido(melhor_global, criterio)
            println("─"^80)
            println("  ✓ CRITÉRIO DE PARADA ATENDIDO na geração $ger")
            break
        end
    end
    println("═"^80, "\n")
    return melhor_global
end

# =========================================================
# EXECUÇÃO E CONFIGURAÇÃO
# =========================================================

drone = Drone(massa = 0.468, gravidade = 9.81, Ixx = 4.856e-3, Iyy = 4.856e-3, Izz = 8.801e-3,
              Ct = 2.980e-6, Cl = 1.14e-7, L = 0.225, estadosIniciais = zeros(12))

config = ConfigAG(40, 100, 0.20, 0.85, 0.15, 0.25, 0.04, 3, 42)
criterio = CriterioParada(10.0, 0.05)
pesos = PesosCusto(2.0, 10.0, 0.5, 2.5)

# Chute Inicial (Modificado para um baseline mais balanceado de matrizes de custo LQR)
genes_iniciais = [
    100.0, 10.0, 1.0,    # Z (Q11, Q22, R)
    100.0, 10.0, 1.0,    # X 
    100.0, 10.0, 1.0,    # Y 
    100.0, 10.0, 1.0,    # Φ (Roll)
    100.0, 10.0, 1.0,    # θ (Pitch)
    100.0, 10.0, 1.0     # Ψ (Yaw)
]

melhor = otimizarLQRAg(drone; config=config, criterio=criterio, pesos=pesos, genes_iniciais=genes_iniciais)

# =========================================================
# RELATÓRIO FINAL
# =========================================================
ctrl_otimo = decodificarGenesLQR(melhor.genes, drone)
params_val = DroneTrajetoriaLQRParams(drone, ctrl_otimo, traj_avaliacao)
sol_final  = resolverSistema(drone_lqr_cascata_mf!, zeros(12), (0.0, T_SIM_AG), params_val; resolucao=0.01)

plotarRespostaFinal(sol_final, criterio)

g = melhor.genes
println("\n╔══════════════════════════════════════════════════════════════════╗")
println("║       PESOS ÓTIMOS (Q E R) ENCONTRADOS PELO AG PARA O LQR        ║")
println("╠══════════════════════════════════════════════════════════════════╣")
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Altitude Z", g[1],  g[2],  g[3])
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Posição X",  g[4],  g[5],  g[6])
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Posição Y",  g[7],  g[8],  g[9])
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Roll Φ",     g[10], g[11], g[12])
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Pitch θ",    g[13], g[14], g[15])
@printf("║  %-10s  Q = diag(%7.1f, %5.1f)    R = [%5.2f]            ║\n", "Yaw Ψ",      g[16], g[17], g[18])
println("╚══════════════════════════════════════════════════════════════════╝\n")

println("# ── Snippet: Atualize os parâmetros de Construção no Script Principal ──")
@printf("ctrl_atitude = construirAtitudeLQR(drone.Ixx, drone.Iyy, drone.Izz;\n")
@printf("    Q_roll  = [%.1f 0.0; 0.0 %.1f], R_roll  = reshape([%.2f], 1, 1),\n", g[10], g[11], g[12])
@printf("    Q_pitch = [%.1f 0.0; 0.0 %.1f], R_pitch = reshape([%.2f], 1, 1),\n", g[13], g[14], g[15])
@printf("    Q_yaw   = [%.1f 0.0; 0.0 %.1f], R_yaw   = reshape([%.2f], 1, 1))\n\n", g[16], g[17], g[18])

@printf("ctrl_altitude = construirAltitudeLQR(drone.massa; gravidade=drone.gravidade,\n")
@printf("    Q_Z = [%.1f 0.0; 0.0 %.1f], R_Z = reshape([%.2f], 1, 1))\n\n", g[1], g[2], g[3])

@printf("ctrl_posicao = construirPosicaoLQR(gravidade=drone.gravidade;\n")
@printf("    Q_X = [%.1f 0.0; 0.0 %.1f], R_X = reshape([%.2f], 1, 1),\n", g[4], g[5], g[6])
@printf("    Q_Y = [%.1f 0.0; 0.0 %.1f], R_Y = reshape([%.2f], 1, 1))\n", g[7], g[8], g[9])