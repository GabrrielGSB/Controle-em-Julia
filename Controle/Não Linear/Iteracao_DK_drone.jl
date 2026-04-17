include("../Funções Auxiliares/SOS.jl")

# =========================================================================
# 1. Funções Auxiliares SOS (Adaptadas para Multi-Input)
# =========================================================================



# =========================================================================
# 2. Definição do Drone Polinomial
# =========================================================================

function construir_drone_polinomial()
    # 12 Estados do Drone
    @polyvar X Y Z u v w_lin phi theta psi p q r
    vars = [X, Y, Z, u, v, w_lin, phi, theta, psi, p, q, r]

    # Parâmetros Físicos
    m = 0.468
    grav = 9.81
    Ixx = 4.856e-3
    Iyy = 4.856e-3
    Izz = 8.801e-3

    # APROXIMAÇÃO DE TAYLOR (Transformando Trig em Polinômios)
    # sin(x) ≈ x ; cos(x) ≈ 1 - x^2/2 ; tan(x) ≈ x
    s_phi = phi;   c_phi = 1 - 0.5*phi^2;  t_theta = theta
    s_theta = theta; c_theta = 1 - 0.5*theta^2
    s_psi = psi;   c_psi = 1 - 0.5*psi^2

    # Construindo vetor f(x) [12x1] - Apenas os termos de até grau 3 mantidos
    f = Vector{Polynomial{true, Float64}}(undef, 12)
    
    f[1] = w_lin * (s_phi * s_psi + c_phi * c_psi * s_theta) - v * (c_phi * s_psi - c_psi * s_phi * s_theta) + u * c_psi * c_theta
    f[2] = v * (c_phi * c_psi + s_phi * s_psi * s_theta) - w_lin * (c_psi * s_phi - c_phi * s_psi * s_theta) + u * c_theta * s_psi
    f[3] = w_lin * c_phi * c_theta - u * s_theta + v * c_theta * s_phi
    
    f[4] = r * v - q * w_lin + grav * s_theta
    f[5] = p * w_lin - r * u - grav * c_theta * s_phi
    f[6] = q * u - p * v - grav * c_phi * c_theta
    
    f[7] = p + q * s_phi * t_theta + r * c_phi * t_theta
    f[8] = q * c_phi - r * s_phi
    f[9] = r * (c_phi / c_theta) + q * (s_phi / c_theta) # Aproximação grossa de divisão
    
    f[10] = ((Iyy - Izz) / Ixx) * q * r
    f[11] = ((Izz - Ixx) / Iyy) * p * r
    f[12] = ((Ixx - Iyy) / Izz) * p * q

    # Construindo matriz g(x) [12x4]
    g = zeros(12, 4)
    g[6, 1]  = 1.0 / m
    g[10, 2] = 1.0 / Ixx
    g[11, 3] = 1.0 / Iyy
    g[12, 4] = 1.0 / Izz

    return f, g, vars
end

# =========================================================================
# 3. Main: Algoritmo D-K para o Drone
# =========================================================================

println("=========================================")
println("Iniciando Iteração D-K para Quadrotor 12-DOF")
println("=========================================\n")

# Extraindo sistema
f, g, vars = construir_drone_polinomial()

# Configurações da Iteração
max_iter = 3
iter = 1

# Chute Inicial V0 (Matriz Identidade ponderada -> V = x^T * Q * x)
# Dando pesos maiores para posição e atitude
pesos = [10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0]
V_atual = sum(pesos[i] * vars[i]^2 for i in 1:12)

u_final = nothing
V_final = V_atual

while iter <= max_iter
    println("--- Iteração $iter ---")
   
    # Passo 1: Achar Controle u(x)
    println("  > Buscando vetor de controle U(x)...")
    u_atual, lamb_atual, success_u = encontrarControleSOS(f, g, V_atual, vars)
    
    if !success_u
        println("  [!] Falha no Passo U. O problema pode ser infactível.")
        break
    end
    
    # Passo 2: Achar Nova Lyapunov V(x)
    println("  > Buscando nova função V(x)...")
    V_novo, success_V = encontrarLyapunovSOS(f, g, u_atual, vars)
    
    if !success_V
        println("  [!] Falha no Passo V. O problema pode ser infactível.")
        break
    end
   
    # Atualiza as variáveis globais
    global V_atual = V_novo 
    global V_final = V_novo
    global u_final = u_atual
    
    println("  [OK] Novo par (U, V) viável encontrado!\n")
    global iter += 1
end

println("=========================================")
if u_final === nothing
    println("O algoritmo falhou em estabilizar o Drone 12-DOF.")
else
    println("Melhor Lei de Controle U(x) encontrada!")
    # O output terá 4 polinômios, representando U1, U2, U3 e U4.
    println("U1 (Empuxo) = ", u_final[1])
    println("U2 (Roll)   = ", u_final[2])
    println("U3 (Pitch)  = ", u_final[3])
    println("U4 (Yaw)    = ", u_final[4])
end