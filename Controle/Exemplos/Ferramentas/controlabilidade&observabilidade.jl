# Controle/Exemplos/Testar_Propriedades.jl

using LinearAlgebra

# =========================================================================
# INCLUDES
    include("../../SistemasLTI/Pêndulo.jl")
    include("../../SistemasLTI/PênduloInvertido.jl") 
    include("../../Ferramentas/LinearizarSistema.jl") 
    include("../../Ferramentas/PropriedadesSistema.jl") 
# =========================================================================

# =========================================================================
# ANÁLISE DO PÊNDULO SIMPLES
    println("==================================================")
    println(" 🛠️  ANÁLISE 1: PÊNDULO SIMPLES")
    println("==================================================")

    pendulo_simples = Pendulo() 

    # x = [ângulo, velocidade_angular]
    x_eq = [0.0, 0.0] 
    u_eq = [0.0]

    h(x, u) = [x[1]]

    A, B, C, D = linearizar(pendulo_simples, x_eq, u_eq, h=h)

    println("Controlabilidade:")
    if controlavel(A, B)
        println("✅ SIM! O motor no eixo consegue controlar o ângulo e a velocidade.")
    else
        println("❌ NÃO!")
    end

    println("\nObservabilidade:")
    if observavel(A, C)
        println("✅ SIM! Medindo apenas o ângulo, o framework deduz a velocidade.")
    else
        println("❌ NÃO!")
    end
# =========================================================================


# =========================================================================
# 3. ANÁLISE DO PÊNDULO INVERTIDO NO CARRINHO
# =========================================================================
println("\n==================================================")
println(" 🛠️  ANÁLISE 2: PÊNDULO INVERTIDO (CART-POLE)")
println("==================================================")

# Instanciando a planta com valores padrão
pendulo_inv = PenduloInvertido() 

# Ponto de equilíbrio instável (haste perfeitamente para cima)
# x = [pos_carrinho, vel_carrinho, angulo_haste, vel_angular_haste]
x_eq = [0.0, 0.0, 0.0, 0.0]
u_eq = [0.0]

# Linearização base sem sensor definido (para pegar A e B)
A, B, _, _ = linearizar(pendulo_inv, x_eq, u_eq)

println("Controlabilidade:")
if controlavel(A, B)
    println("   ✅ SIM! A força no carrinho é suficiente para estabilizar tudo.")
else
    println("   ❌ NÃO.")
end

# --- TESTE DE SENSORES (OBSERVABILIDADE) ---

println("\nObservabilidade: Testando diferentes sensores...")

# Cenário A: Medindo APENAS a Posição do Carrinho (Estado 1)
h_carrinho(x, u) = [x[1]]
_, _, C_car, _ = linearizar(pendulo_inv, x_eq, u_eq, h=h_carrinho)

print("   Sensor no Carrinho: ")
if observavel(A_inv, C_car)
    println("✅ OBSERVÁVEL (O movimento do pêndulo afeta o carrinho, logo deduzimos tudo).")
else
    println("❌ NÃO OBSERVÁVEL.")
end

# Cenário B: Medindo APENAS o Ângulo da Haste (Estado 3)
h_angulo(x, u) = [x[3]]
_, _, C_ang, _ = linearizar(pendulo_inv, x_eq, u_eq, h=h_angulo)

print("   Sensor na Haste:    ")
if observavel(A_inv, C_ang)
    println("✅ OBSERVÁVEL.")
else
    println("❌ NÃO OBSERVÁVEL (A posição absoluta X do carrinho não altera o ângulo do pêndulo. O sistema perde a noção de onde o carrinho está na pista!).")
end