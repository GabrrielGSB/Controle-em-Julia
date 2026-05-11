include("../../SistemasLTI/PênduloInvertido.jl")

# 1. Instanciamos a planta física do seu framework para carregar os parâmetros (M, m, l, g, b, etc.)
pendulo = PenduloInvertido()

x_eq = [0.0, 0.0, 0.0, 0.0]
u_eq = [0.0]

A, B, C, D = linearizar_sistema(pendulo, x_eq, u_eq)

println("=== Matrizes do Espaço de Estados Linearizado ===")
println("\nMatriz A (Dinâmica dos Estados):")
display(A)

println("\nMatriz B (Matriz de Entrada):")
display(B)

println("\nMatriz C (Matriz de Saída):")
display(C)

polos = eigvals(A)

println("Os polos do sistema (Autovalores de A) são:")
display(polos)