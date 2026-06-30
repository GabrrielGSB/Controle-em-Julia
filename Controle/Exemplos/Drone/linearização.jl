# ===================================================================
# linearizar_drone.jl
# 
# Exemplo de linearização do drone quadrotor 12-DOF no ponto de 
# equilíbrio (pairando). Obtém as matrizes A, B, C, D do sistema 
# linearizado e as exibe de forma organizada no terminal.
# ===================================================================

# 1. INCLUDES NECESSÁRIOS
include("../../SistemasLTI/Drone.jl")
include("../../Ferramentas/LinearizarSistema.jl")

# 2. PARÂMETROS FÍSICOS DO DRONE
m   = 0.468      # massa [kg]
g   = 9.81       # gravidade [m/s²]
Ixx = 4.856e-3
Iyy = 4.856e-3
Izz = 8.801e-3
Ct  = 2.980e-6
Cl  = 1.14e-7
L   = 0.225

# Instancia a planta (estados iniciais não importam para a linearização,
# pois o ponto de operação será sobrescrito)
drone = Drone(
    massa       = m,
    gravidade   = g,
    Ixx         = Ixx,
    Iyy         = Iyy,
    Izz         = Izz,
    Ct          = Ct,
    Cl          = Cl,
    L           = L,
    estadosIniciais = zeros(12)
)

# 3. PONTO DE EQUILÍBRIO (VOO PAIRADO)
# Estados: posição (0,0,0), velocidades lineares nulas, ângulos nulos, 
#          velocidades angulares nulas.
x0 = zeros(12)

# Entradas: para pairar, U1 deve igualar o peso (m*g). Rols, pitch e yaw = 0.
u0 = [m * g, 0.0, 0.0, 0.0]  # [U1, U2, U3, U4]

println("\n" * "="^70)
println("  LINEARIZAÇÃO DO DRONE NO PONTO DE EQUILÍBRIO (HOVER)")
println("  x0 = zeros(12)  |  u0 = [$(round(u0[1], digits=2)), 0, 0, 0]")
println("="^70)

# ===================================================================
# 4. CASO 1: SAÍDA COMPLETA (C = I, D = 0)
#    Quando não passamos a função h(x,u), a linearização assume que
#    medimos todos os estados.
# ===================================================================
A, B, C, D = linearizar(drone, x0, u0)

println("\n▶ MATRIZES DO SISTEMA COM SAÍDA COMPLETA (y = x)")

# 4.1 Função auxiliar para impressão bonita
function imprimir_matriz(M, nome; precisao=4)
    println("\n" * "-"^60)
    println("📌 $nome  |  Dimensão: $(size(M,1)) × $(size(M,2))")
    println("-"^60)
    M_arredondada = round.(M, digits=precisao)
    display(M_arredondada)
end

imprimir_matriz(A, "Matriz A (Dinâmica dos Estados)")
imprimir_matriz(B, "Matriz B (Entrada de Controle)")
imprimir_matriz(C, "Matriz C (Saída)")
imprimir_matriz(D, "Matriz D (Alimentação direta)")

# ===================================================================
# 5. CASO 2: SAÍDA PERSONALIZADA (MEDINDO APENAS Z, Φ, θ, Ψ)
#    Útil para projetar controladores que usam apenas sensores de 
#    altitude (barômetro) e atitude (IMU).
# ===================================================================
# Definição da função de saída: y = h(x, u) = [Z, Φ, θ, Ψ]
function h_drone(x, u)
    # x[3] = Z, x[7] = Φ (Roll), x[8] = θ (Pitch), x[9] = Ψ (Yaw)
    return [x[3], x[7], x[8], x[9]]
end

A_custom, B_custom, C_custom, D_custom = linearizar(drone, x0, u0; h=h_drone)

println("\n\n▶ MATRIZES COM SAÍDA PERSONALIZADA: y = [Z, Φ, θ, Ψ]")

imprimir_matriz(A_custom, "Matriz A (inalterada)")
imprimir_matriz(B_custom, "Matriz B (inalterada)")
imprimir_matriz(C_custom, "Matriz C (Saída customizada)")
imprimir_matriz(D_custom, "Matriz D (Saída customizada)")

# ===================================================================
# 6. ANÁLISE RÁPIDA DOS POLOS (ESTABILIDADE DO LINEARIZADO)
# ===================================================================
using LinearAlgebra

polos = eigvals(A)
println("\n" * "="^70)
println("  POLOS DO SISTEMA LINEARIZADO (A)")
println("-"^70)
# Exibe apenas a parte real e imaginária para ficar mais limpo
for (i, p) in enumerate(polos)
    println("  Polo $(lpad(i,2)):  $(round(real(p), digits=4))  $(round(imag(p), digits=4))im")
end

# Contagem de polos instáveis (parte real > 0)
instaveis = count(real(polos) .> 1e-6)
println("-"^70)
if instaveis > 0
    println("  ⚠️  O sistema linearizado possui $instaveis polo(s) instável(is).")
    println("     (Esperado, pois o drone em hover é instável em atitude e posição).")
else
    println("  ✅ Sistema linearizado estável (todos os polos com Re < 0).")
end
println("="^70)