using Symbolics, LinearAlgebra, ControlSystems

# ========== 1. Definição simbólica do sistema ==========
@variables t
# Estados (12)
x, y, z, φ, θ, ψ, ωx, ωy, ωz, vx, vy, vz = states = [x(t), y(t), z(t), φ(t), θ(t), ψ(t), ωx(t), ωy(t), ωz(t), vx(t), vy(t), vz(t)]
# Entradas (4): torques e empuxo total
@variables τφ(t) τθ(t) τψ(t) T(t)
# Parâmetros
@variables m g Jxx Jyy Jzz

# Funções trigonométricas
cφ = cos(φ); sφ = sin(φ); cθ = cos(θ); sθ = sin(θ); cψ = cos(ψ); sψ = sin(ψ); tθ = tan(θ)

# Derivadas (mesmo do seu arquivo)
dx = (cθ*cψ)*vx + (sφ*sθ*cψ - cφ*sψ)*vy + (cφ*sθ*cψ + sφ*sψ)*vz
dy = (cθ*sψ)*vx + (sφ*sθ*sψ + cφ*cψ)*vy + (cφ*sθ*sψ - sφ*cψ)*vz
dz = (-sθ)*vx + (sφ*cθ)*vy + (cφ*cθ)*vz

dφ = ωx + sφ*tθ*ωy + cφ*tθ*ωz
dθ = cφ*ωy - sφ*ωz
dψ = (sφ/cθ)*ωy + (cφ/cθ)*ωz

dωx = (1/Jxx)*(τφ - (Jzz - Jyy)*ωy*ωz)
dωy = (1/Jyy)*(τθ - (Jxx - Jzz)*ωz*ωx)
dωz = (1/Jzz)*(τψ - (Jyy - Jxx)*ωx*ωy)

dvx = ωz*vy - ωy*vz + g*sθ
dvy = ωz*vx + ωx*vz - g*cθ*sφ
dvz = ωy*vx - ωx*vy - g*cφ*cθ + T/m

f = [dx, dy, dz, dφ, dθ, dψ, dωx, dωy, dωz, dvx, dvy, dvz]

# ========== 2. Linearização no hover ==========
X = [x, y, z, φ, θ, ψ, ωx, ωy, ωz, vx, vy, vz]
U = [τφ, τθ, τψ, T]

A_sym = Symbolics.jacobian(f, X)        # 12x12
B_sym = Symbolics.jacobian(f, U)        # 12x4

# Parâmetros numéricos (exemplo)
p = Dict(m => 1.0, g => 9.81, Jxx => 0.01, Jyy => 0.01, Jzz => 0.02)

# Ponto de equilíbrio: estados zero, entradas constantes
eq = Dict(
    x=>0, y=>0, z=>0, φ=>0, θ=>0, ψ=>0,
    ωx=>0, ωy=>0, ωz=>0, vx=>0, vy=>0, vz=>0,
    τφ=>0, τθ=>0, τψ=>0, T=>p[m]*p[g]
)

# Substituir e converter para matrizes numéricas
A = float.(Symbolics.value.(Symbolics.substitute.(A_sym, (eq, p))))
B = float.(Symbolics.value.(Symbolics.substitute.(B_sym, (eq, p))))

println("Matriz A (sem controle):")
display(A)

# ========== 3. Projetar controlador LQR ==========
# Matrizes de peso (ajuste fino)
Q_lqr = diagm(ones(12))        # peso nos estados
R_lqr = 10.0 * I(4)            # peso nas entradas (reduz esforço)

# Resolver ganho LQR (sistema linear contínuo)
K = lqr(A, B, Q_lqr, R_lqr)
Acl = A - B * K

println("\nMatriz de malha fechada Acl = A - BK:")
display(Acl)

# ========== 4. Encontrar P pela equação de Lyapunov ==========
Q = I(12)                       # definida positiva
P = lyap(Acl', Q)               # resolve Acl' * P + P * Acl = -Q

# Verificar se P é definida positiva (todos os autovalores > 0)
eigP = eigvals(P)
println("\nAutovalores de P:")
println(eigP)

if all(eigP .> 0)
    println("\n✅ P é definida positiva → V(x) = xᵀP x é função de Lyapunov.")
    println("✅ Derivada: V̇(x) = -xᵀQ x < 0 para x≠0.")
    println("→ O sistema em malha fechada é localmente assintoticamente estável.")
else
    println("\n❌ P não é definida positiva → a função candidata não prova estabilidade.")
end

# ========== 5. (Opcional) Verificar instabilidade sem controle ==========
# Tentar resolver para A (sem controle)
P_instab = lyap(A', Q)   # pode não ser definida positiva
eigP_instab = eigvals(P_instab)
println("\nAutovalores de P para o sistema sem controle:")
println(eigP_instab)

if any(eigP_instab .< 0)
    println("❌ Existem autovalores negativos → não é possível provar estabilidade com esta Q.")
    println("De fato, o sistema sem controle é instável (autovalores de A com parte real positiva).")
end