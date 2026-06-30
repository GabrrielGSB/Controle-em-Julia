# Controle/Exemplos/Drone/Drone_MPC_Cascata_Trajetoria.jl
#
# SEGUIMENTO DE TRAJETÓRIA COM MPC EM CASCATA — QUADROTOR
# (Utilizando Solver Profissional OSQP.jl)

# =========================================================
# INCLUDES
    include("../../SistemasLTI/Drone.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

using PlotlyJS
using Printf
using LinearAlgebra
using Statistics
using OSQP          # Solver Profissional de QP
using SparseArrays  # Exigência do OSQP

# =========================================================
# TRAJETÓRIAS DE REFERÊNCIA 
# =========================================================

function trajetoriaHover(t; altitude=2.0, t_subida=5.0)
    τ = clamp(t / t_subida, 0.0, 1.0)
    h = altitude * τ^2 * (3 - 2τ)
    return (X=0.0, Y=0.0, Z=h, Ψ=0.0)
end

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

# =========================================================
# PARÂMETROS DO MPC (PGDM Removido)
# =========================================================

struct ParamsMPC_Pos
    Np     ::Int
    Nc     ::Int
    dt     ::Float64
    Q      ::Matrix{Float64}
    R      ::Matrix{Float64}
    P      ::Matrix{Float64}
    u_min  ::Vector{Float64}
    u_max  ::Vector{Float64}
    Δu_max ::Vector{Float64}
end

struct ParamsMPC_Att
    Np     ::Int
    Nc     ::Int
    dt     ::Float64
    Q      ::Matrix{Float64}
    R      ::Matrix{Float64}
    P      ::Matrix{Float64}
    u_min  ::Vector{Float64}
    u_max  ::Vector{Float64}
    Δu_max ::Vector{Float64}
end

struct DroneMPCParams
    drone      ::Drone
    mpc_pos    ::ParamsMPC_Pos
    mpc_att    ::ParamsMPC_Att
    trajetoria ::Function
end

# =========================================================
# LINEARIZAÇÃO E CONSTRUÇÃO MATRICIAL DO MODELO
# =========================================================

function linearizar_pos(x::AbstractVector, u_op::AbstractVector, d::Drone, dt::Float64)
    g = d.gravidade; m = d.massa; Ψ = x[9]

    Ac = zeros(6, 6)
    Ac[1, 4] = 1.0; Ac[2, 5] = 1.0; Ac[3, 6] = 1.0

    Bc = zeros(6, 3)
    Bc[4, 2] =  g * cos(Ψ); Bc[4, 3] =  g * sin(Ψ)
    Bc[5, 2] =  g * sin(Ψ); Bc[5, 3] = -g * cos(Ψ)
    Bc[6, 1] =  1.0 / m

    Ad = I(6) + Ac * dt
    Bd = Bc * dt
    return Ad, Bd
end

function linearizar_att(x::AbstractVector, u_op::AbstractVector, d::Drone, dt::Float64)
    Ixx = d.Ixx; Iyy = d.Iyy; Izz = d.Izz

    Ac = zeros(6, 6)
    Ac[1, 4] = 1.0; Ac[2, 5] = 1.0; Ac[3, 6] = 1.0

    Bc = zeros(6, 3)
    Bc[4, 1] = 1.0 / Ixx; Bc[5, 2] = 1.0 / Iyy; Bc[6, 3] = 1.0 / Izz

    Ad = I(6) + Ac * dt
    Bd = Bc * dt
    return Ad, Bd
end

function construir_matrizes_mpc(Ad::Matrix, Bd::Matrix, Q::Matrix, R::Matrix, P::Matrix, Np::Int, Nc::Int)
    nx = size(Ad, 1); nu = size(Bd, 2)

    Φ_mat = zeros(Np * nx, nx)
    Ak    = Matrix{Float64}(I, nx, nx)
    for k = 1:Np
        Ak = Ad * Ak
        Φ_mat[(k-1)*nx+1 : k*nx, :] = Ak
    end

    Γ_mat = zeros(Np * nx, Nc * nu)
    for i = 1:Np
        for j = 1:min(i, Nc)
            pot  = i - j
            Apot = Matrix{Float64}(I, nx, nx)
            for _ = 1:pot; Apot = Ad * Apot; end
            Γ_mat[(i-1)*nx+1 : i*nx, (j-1)*nu+1 : j*nu] = Apot * Bd
        end
    end

    Q_bar = zeros(Np * nx, Np * nx)
    for k = 1:Np
        bloco = (k < Np) ? Q : P
        Q_bar[(k-1)*nx+1 : k*nx, (k-1)*nx+1 : k*nx] = bloco
    end
    R_bar = kron(I(Nc), R)

    H = Γ_mat' * Q_bar * Γ_mat + R_bar
    H = (H + H') / 2.0  # Simetria numérica estrita

    return Φ_mat, Γ_mat, H, Q_bar, R_bar
end

# =========================================================
# NOVO OTIMIZADOR: OSQP SOLVER (Substituindo o PGDM)
# =========================================================

"""
    osqp_solve(H, f, ΔU_lb, ΔU_ub)

Resolve o problema QP usando a biblioteca profissional OSQP.
Minimiza: (1/2) x' H x + f' x
Sujeito a: lb ≤ x ≤ ub
"""
function osqp_solve(H::Matrix, f::Vector, lb::Vector, ub::Vector)
    n = length(f)
    
    # OSQP exige que as matrizes passem pelo formato sparse do Julia
    P_sparse = sparse(H)
    q_vec    = f
    
    # A restrição é sobre a própria variável de decisão: I * ΔU
    A_sparse = sparse(I, n, n)
    
    # Instancia e configura o modelo silenciosamente
    model = OSQP.Model()
    OSQP.setup!(model; P=P_sparse, q=q_vec, A=A_sparse, l=lb, u=ub, 
                verbose=false, eps_abs=1e-4, eps_rel=1e-4)
    
    # Resolve
    results = OSQP.solve!(model)
    
    # Se falhar em convergir por algum motivo numérico grave, retorna zero 
    # para evitar que o solver numérico EDO receba lixo e quebre a simulação.
    if results.info.status_val != 1 && results.info.status_val != 2
        return zeros(n)
    end
    
    return results.x
end

# =========================================================
# MPC DE POSIÇÃO (MALHA EXTERNA)
# =========================================================

function mpc_pos_step(x_drone::Vector, u_prev_pos::Vector, ref_traj::Function, t::Float64, p::DroneMPCParams)
    d = p.drone; mp = p.mpc_pos; dt = mp.dt; Np = mp.Np; Nc = mp.Nc
    g = d.gravidade; m = d.massa

    Ψ = x_drone[9]; Vx = x_drone[4]; Vy = x_drone[5]; Vz = x_drone[6]
    Vx_i =  Vx * cos(Ψ) - Vy * sin(Ψ)
    Vy_i =  Vx * sin(Ψ) + Vy * cos(Ψ)

    z0 = [x_drone[1], x_drone[2], x_drone[3], Vx_i, Vy_i, Vz]

    Ad, Bd = linearizar_pos(x_drone, u_prev_pos, d, dt)

    nx = 6; nu = 3
    X_ref = zeros(Np * nx)
    for k = 1:Np
        ref  = ref_traj(t + k * dt)
        dref = ref_traj(t + (k+1) * dt)
        Vx_r = (dref.X - ref.X) / dt
        Vy_r = (dref.Y - ref.Y) / dt
        Vz_r = (dref.Z - ref.Z) / dt
        X_ref[(k-1)*nx+1 : k*nx] = [ref.X, ref.Y, ref.Z, Vx_r, Vy_r, Vz_r]
    end

    Φ_mat, Γ_mat, H, Q_bar, _ = construir_matrizes_mpc(Ad, Bd, mp.Q, mp.R, mp.P, Np, Nc)

    livre = Φ_mat * z0
    erro  = livre - X_ref
    f_vec = Γ_mat' * Q_bar * erro

    ΔU_lb = repeat(max.(mp.u_min .- u_prev_pos, -mp.Δu_max), Nc)
    ΔU_ub = repeat(min.(mp.u_max .- u_prev_pos,  mp.Δu_max), Nc)

    # Chamada do Solver Profissional
    ΔU_sol = osqp_solve(H, f_vec, ΔU_lb, ΔU_ub)

    Δu1 = ΔU_sol[1:nu]
    u_opt = clamp.(u_prev_pos + Δu1, mp.u_min, mp.u_max)

    return u_opt[1], u_opt[2], u_opt[3], u_opt
end

# =========================================================
# MPC DE ATITUDE (MALHA INTERNA)
# =========================================================

function mpc_att_step(x_drone::Vector, u_prev_att::Vector, Φ_des::Float64, θ_des::Float64, Ψ_r::Float64, t::Float64, p::DroneMPCParams)
    d = p.drone; ma = p.mpc_att; dt = ma.dt; Np = ma.Np; Nc = ma.Nc

    Φ = x_drone[7]; θ = x_drone[8]; Ψ = x_drone[9]
    VΦ = x_drone[10]; Vθ = x_drone[11]; VΨ = x_drone[12]

    e_Ψ_0 = atan(sin(Ψ_r - Ψ), cos(Ψ_r - Ψ))
    z0    = [Φ, θ, e_Ψ_0 + Ψ, VΦ, Vθ, VΨ]

    Ad, Bd = linearizar_att(x_drone, u_prev_att, d, dt)

    nx = 6; nu = 3
    X_ref = zeros(Np * nx)
    for k = 1:Np
        X_ref[(k-1)*nx+1 : k*nx] = [Φ_des, θ_des, Ψ_r, 0.0, 0.0, 0.0]
    end

    Φ_mat, Γ_mat, H, Q_bar, _ = construir_matrizes_mpc(Ad, Bd, ma.Q, ma.R, ma.P, Np, Nc)

    livre = Φ_mat * z0
    erro  = livre - X_ref
    f_vec = Γ_mat' * Q_bar * erro

    ΔU_lb = repeat(max.(ma.u_min .- u_prev_att, -ma.Δu_max), Nc)
    ΔU_ub = repeat(min.(ma.u_max .- u_prev_att,  ma.Δu_max), Nc)

    # Chamada do Solver Profissional
    ΔU_sol = osqp_solve(H, f_vec, ΔU_lb, ΔU_ub)

    Δu1 = ΔU_sol[1:nu]
    u_opt = clamp.(u_prev_att + Δu1, ma.u_min, ma.u_max)

    return u_opt[1], u_opt[2], u_opt[3], u_opt
end

# =========================================================
# DINÂMICA DA MALHA FECHADA MPC EM CASCATA
# =========================================================

mutable struct DroneMPCState
    u_pos ::Vector{Float64}
    u_att ::Vector{Float64}
end

function drone_mpc_cascata_mf!(dx, x, p::Tuple{DroneMPCParams, DroneMPCState}, t)
    params, mpc_state = p
    d = params.drone
    g = d.gravidade
    m = d.massa

    x_drone = x[1:12]

    ref = params.trajetoria(t)
    Ψ_r = ref.Ψ

    U1, θ_des, Φ_des, u_pos_new = mpc_pos_step(x_drone, mpc_state.u_pos, params.trajetoria, t, params)
    U1_total = clamp(U1 + m * g, 0.0, 2.5 * m * g)

    U2, U3, U4, u_att_new = mpc_att_step(x_drone, mpc_state.u_att, Φ_des, θ_des, Ψ_r, t, params)

    mpc_state.u_pos = u_pos_new
    mpc_state.u_att = u_att_new

    d.dinamica!(view(dx, 1:12), x_drone, d, t; u=[U1_total, U2, U3, U4])
end

# =========================================================
# EXECUÇÃO E CONFIGURAÇÃO
# =========================================================

drone = Drone(massa = 0.468, gravidade = 9.81, Ixx = 4.856e-3, Iyy = 4.856e-3, Izz = 8.801e-3,
              Ct = 2.980e-6, Cl = 1.14e-7, L = 0.225, estadosIniciais = zeros(12))

Q_pos = Diagonal([5.0, 5.0, 8.0, 1.0, 1.0, 2.0])
R_pos = Diagonal([0.01, 5.0, 5.0])
P_pos = 2.0 * Q_pos
mpc_pos = ParamsMPC_Pos(10, 4, 0.05, Matrix(Q_pos), Matrix(R_pos), Matrix(P_pos),
                        [0.0, -0.35, -0.35], [20.0, 0.35, 0.35], [2.0, 0.08, 0.08])

Q_att = Diagonal([20.0, 20.0, 15.0, 2.0, 2.0, 1.5])
R_att = Diagonal([0.5, 0.5, 0.5])
P_att = 3.0 * Q_att
mpc_att = ParamsMPC_Att(6, 3, 0.02, Matrix(Q_att), Matrix(R_att), Matrix(P_att),
                        [-0.5, -0.5, -0.2], [0.5, 0.5, 0.2], [0.15, 0.15, 0.08])

tspan = (0.0, 40.0)
trajetoria = t -> trajetoriaHelice(t; raio=1.5, ω_giro=0.4, v_subida=0.08, t_espera=4.0)
nome_traj = "Hélice Ascendente (OSQP)"

println("="^60)
println("  DRONE — MPC EM CASCATA (Com OSQP Solver)")
println("  Aviso: Simulação pode demorar mais devido ao overhead do OSQP no EDO.")
println("="^60)

params    = DroneMPCParams(drone, mpc_pos, mpc_att, trajetoria)
mpc_state = DroneMPCState([drone.massa * drone.gravidade, 0.0, 0.0], [0.0, 0.0, 0.0])

mf! = (dx, x, p, t) -> drone_mpc_cascata_mf!(dx, x, p, t)
x0 = zeros(12)

# Tempo de resolução aqui será superior, mas com matemática estritamente correta.
sol = resolverSistema(mf!, x0, tspan, (params, mpc_state); resolucao = 0.02)


function plotarTrajetoria3D(solucao, trajetoria, tspan; titulo="Seguimento de Trajetória 3D")
    t_vec  = solucao.t
    X_real = [u[1] for u in solucao.u]
    Y_real = [u[2] for u in solucao.u]
    Z_real = [u[3] for u in solucao.u]

    refs  = [trajetoria(t) for t in t_vec]
    X_ref = [r.X for r in refs]
    Y_ref = [r.Y for r in refs]
    Z_ref = [r.Z for r in refs]

    trace_real  = scatter3d(x=X_real, y=Y_real, z=Z_real, mode="lines",
                             name="Drone (MPC)", line=attr(color="#1f77b4", width=4))
    trace_ref   = scatter3d(x=X_ref,  y=Y_ref,  z=Z_ref,  mode="lines",
                             name="Referência",  line=attr(color="#d62728", width=3, dash="dash"))
    trace_inicio = scatter3d(x=[X_real[1]], y=[Y_real[1]], z=[Z_real[1]],
                              mode="markers", name="Início", marker=attr(size=6, color="green"))
    trace_fim    = scatter3d(x=[X_real[end]], y=[Y_real[end]], z=[Z_real[end]],
                              mode="markers", name="Fim",
                              marker=attr(size=6, color="red", symbol="x"))

    layout = Layout(title_text=titulo, title_x=0.5,
                    scene=attr(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Z (m)"),
                    width=850, height=650)
    display(plot([trace_real, trace_ref, trace_inicio, trace_fim], layout))
end

function plotarErrosRastreamento(solucao, trajetoria; titulo="Erros de Rastreamento — MPC")
    t_vec  = solucao.t
    refs   = [trajetoria(t) for t in t_vec]
    e_X    = [refs[i].X - solucao.u[i][1] for i in eachindex(t_vec)]
    e_Y    = [refs[i].Y - solucao.u[i][2] for i in eachindex(t_vec)]
    e_Z    = [refs[i].Z - solucao.u[i][3] for i in eachindex(t_vec)]
    erro_3d = [sqrt(e_X[i]^2 + e_Y[i]^2 + e_Z[i]^2) for i in eachindex(t_vec)]

    fig = make_subplots(rows=4, cols=1, shared_xaxes=true, vertical_spacing=0.08,
                        subplot_titles=reshape(["Erro X (m)", "Erro Y (m)",
                                                "Erro Z (m)", "Erro Euclidiano 3D (m)"], 1, 4))
    for (linha, (dados, cor, nome)) in enumerate([
            (e_X,    "#1f77b4", "eₓ"),
            (e_Y,    "#2ca02c", "e_Y"),
            (e_Z,    "#d62728", "e_Z"),
            (erro_3d,"#9467bd", "‖e‖₂")])
        add_trace!(fig, scatter(x=t_vec, y=dados, mode="lines",
                                name=nome, line=attr(color=cor, width=2)), row=linha, col=1)
        add_trace!(fig, scatter(x=[t_vec[1], t_vec[end]], y=[0.0, 0.0],
                                mode="lines", showlegend=false,
                                line=attr(color="black", width=1, dash="dot")), row=linha, col=1)
    end
    relayout!(fig, title_text=titulo, title_x=0.5,
              height=800, width=800, hovermode="x unified")
    relayout!(fig, xaxis4_title="Tempo (s)")
    display(fig)
end

function plotarAtitudeEmpuxo(solucao)
    t_vec   = solucao.t
    Φ_hist  = [u[7] for u in solucao.u]
    θ_hist  = [u[8] for u in solucao.u]
    Ψ_hist  = [u[9] for u in solucao.u]
    Vz_hist = [u[6] for u in solucao.u]
    Z_hist  = [u[3] for u in solucao.u]

    fig = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.1,
                        subplot_titles=reshape(["Altitude Z (m) e Velocidade Vz (m/s)",
                                                "Roll Φ e Pitch θ (rad)",
                                                "Yaw Ψ (rad)"], 1, 3))

    add_trace!(fig, scatter(x=t_vec, y=Z_hist,  name="Z",
                            line=attr(color="#1f77b4", width=2)), row=1, col=1)
    add_trace!(fig, scatter(x=t_vec, y=Vz_hist, name="Vz",
                            line=attr(color="#aec7e8", width=2, dash="dot")), row=1, col=1)
    add_trace!(fig, scatter(x=t_vec, y=Φ_hist,  name="Φ (Roll)",
                            line=attr(color="#d62728", width=2)), row=2, col=1)
    add_trace!(fig, scatter(x=t_vec, y=θ_hist,  name="θ (Pitch)",
                            line=attr(color="#ff7f0e", width=2)), row=2, col=1)
    add_trace!(fig, scatter(x=t_vec, y=Ψ_hist,  name="Ψ (Yaw)",
                            line=attr(color="#2ca02c", width=2)), row=3, col=1)

    relayout!(fig, title_text="Atitude e Altitude — Drone MPC Cascata",
              title_x=0.5, height=700, width=800, hovermode="x unified")
    relayout!(fig, xaxis3_title="Tempo (s)")
    display(fig)
end


# ── Métricas de rastreamento ─────────────────────────────────────
refs_final = [trajetoria(t) for t in sol.t]
e_X_final  = [refs_final[i].X - sol.u[i][1] for i in eachindex(sol.t)]
e_Y_final  = [refs_final[i].Y - sol.u[i][2] for i in eachindex(sol.t)]
e_Z_final  = [refs_final[i].Z - sol.u[i][3] for i in eachindex(sol.t)]
rmse_3d    = sqrt(mean(e_X_final.^2 .+ e_Y_final.^2 .+ e_Z_final.^2))

@printf("\n  RMSE de rastreamento 3D : %.4f m\n",   rmse_3d)
@printf("  Erro final em X         : %.4f m\n",   e_X_final[end])
@printf("  Erro final em Y         : %.4f m\n",   e_Y_final[end])
@printf("  Erro final em Z         : %.4f m\n\n", e_Z_final[end])

Φ_max = maximum(abs(sol.u[i][7]) for i in eachindex(sol.t))
θ_max = maximum(abs(sol.u[i][8]) for i in eachindex(sol.t))
@printf("  |Φ|_max (Roll)  : %.2f°\n", rad2deg(Φ_max))
@printf("  |θ|_max (Pitch) : %.2f°\n", rad2deg(θ_max))
println("="^60, "\n")

# =========================================================
# VISUALIZAÇÃO
# =========================================================

plotarTrajetoria3D(sol, trajetoria, tspan;
                   titulo = "Drone MPC Cascata — $nome_traj")

plotarErrosRastreamento(sol, trajetoria;
                        titulo = "Erros de Rastreamento MPC — $nome_traj")

plotarAtitudeEmpuxo(sol)

# =========================================================
# ANÁLISE DE PERFORMANCE
# =========================================================

m_Z = analisarPerformance(sol;
    referencia = trajetoria(tspan[2]).Z,
    idx_estado = 3,
    banda      = 0.05)

imprimirRelatorio(m_Z, nome = "Altitude Z — MPC Cascata")