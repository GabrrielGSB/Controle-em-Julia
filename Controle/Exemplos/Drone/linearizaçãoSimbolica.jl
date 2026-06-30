# =====================================================================
# linearizar_drone_simbolico.jl
# Linearização simbólica do drone 12-DOF usando DynamicPolynomials.
# =====================================================================

using DynamicPolynomials
using LinearAlgebra

# 1. OBTÉM O MODELO POLINOMIAL DO DRONE
include("../../SistemasLTI/Drone.jl")

f, g, vars = obterEqPolyDrone()

# -- DEPURAÇÃO: verifique os tipos retornados --
println("Tipo de f: ", typeof(f))
println("Tipo de g: ", typeof(g))
println("Tipo de vars: ", typeof(vars))
println("Tipo de f[1]: ", typeof(f[1]))
println("Tipo de g[1,1]: ", typeof(g[1,1]))

# 2. DEFINE AS VARIÁVEIS DE ENTRADA SIMBÓLICAS
@polyvar U1 U2 U3 U4
u_vec = [U1, U2, U3, U4]

# 3. CONVERTE g PARA MATRIZ DE POLINÔMIOS CONSTANTES
#    Isso garante que a multiplicação g*u retorne polinômios.
g_poly = [g[i,j] * one(U1) for i in 1:12, j in 1:4]

# 4. CONSTRÓI A DINÂMICA COMPLETA: dx = f(x) + g(x)*u
dx = Vector{Polynomial{true, Float64}}(undef, 12)
for i in 1:12
    # Soma sobre as entradas: ∑ g[i,j] * u_j
    soma = sum(g_poly[i,j] * u_vec[j] for j in 1:4)
    dx[i] = f[i] + soma
end

# 5. CALCULA AS MATRIZES JACOBIANAS SIMBÓLICAS
A_sym = [differentiate(dx[i], vars[j]) for i in 1:12, j in 1:12]
B_sym = [differentiate(dx[i], u_vec[j]) for i in 1:12, j in 1:4]

# 6. DEFINE A SAÍDA (C e D) — Medindo apenas Z, Φ, θ, Ψ
h_sym = [vars[3], vars[7], vars[8], vars[9]]   # [Z, phi, theta, psi]
C_sym = [differentiate(h_sym[i], vars[j]) for i in 1:4, j in 1:12]
D_sym = zeros(4, 4)   # saída não depende diretamente da entrada

# =====================================================================
# 7. IMPRESSÃO DAS MATRIZES SIMBÓLICAS (ENTRADAS NÃO-NULAS)
# =====================================================================
function print_nonzero(M, nome)
    println("\n📌 $nome (Simbólico) — Entradas não-nulas:")
    for i in 1:size(M,1), j in 1:size(M,2)
        expr = M[i,j]
        # Verifica se a expressão não é zero (polinômio nulo)
        if !iszero(expr)
            println("    M[$i,$j] = $expr")
        end
    end
end

println("=" * 70)
println("  LINEARIZAÇÃO SIMBÓLICA DO DRONE")
println("  Modelo polinomial (Taylor) extraído de obterEqPolyDrone()")
println("=" * 70)

print_nonzero(A_sym, "Matriz A (∂f/∂x)")
print_nonzero(B_sym, "Matriz B (∂f/∂u)")
print_nonzero(C_sym, "Matriz C (∂h/∂x)")

# =====================================================================
# 8. SUBSTITUIÇÃO NO PONTO DE EQUILÍBRIO (HOVER)
#    x0 = 0,  u0 = [m*g, 0, 0, 0]
# =====================================================================
m = 0.468
g0 = 9.81

# Dicionário único com todas as variáveis (estados + entradas)
todas_var = [vars; u_vec]
valores = [zeros(12); [m*g0, 0.0, 0.0, 0.0]]
dict_sub = Dict(todas_var[i] => valores[i] for i in 1:16)

# Para A e C (dependem só de x), usamos apenas as substituições das variáveis de estado
dict_x = Dict(vars[i] => 0.0 for i in 1:12)
dict_u = Dict(u_vec[j] => [m*g0,0,0,0][j] for j in 1:4)

# Combinar os dicionários (as entradas não aparecem em A e C, mas isso é seguro)
dict_completo = merge(dict_x, dict_u)

A_num = Float64.(substitute.(A_sym, (dict_completo,)))
C_num = Float64.(substitute.(C_sym, (dict_completo,)))
B_num = Float64.(substitute.(B_sym, (dict_completo,)))
D_num = D_sym   # já é zero

# =====================================================================
# 9. IMPRESSÃO DAS MATRIZES NUMÉRICAS (FORMATAÇÃO BONITA)
# =====================================================================
function imprimir_matriz(M, nome; precisao=4)
    println("\n" * "-" * 60)
    println("📌 $nome  |  Dimensão: $(size(M,1)) × $(size(M,2))")
    println("-" * 60)
    M_arred = round.(M, digits=precisao)
    display(M_arred)
end

println("\n\n▶ MATRIZES AVALIADAS NO PONTO DE EQUILÍBRIO (HOVER)")
imprimir_matriz(A_num, "Matriz A (numérica)")
imprimir_matriz(B_num, "Matriz B (numérica)")
imprimir_matriz(C_num, "Matriz C (numérica)")
imprimir_matriz(D_num, "Matriz D (numérica)")

# =====================================================================
# 10. POLOS DO SISTEMA LINEARIZADO
# =====================================================================
polos = eigvals(A_num)

println("\n" * "=" * 70)
println("  POLOS DO SISTEMA LINEARIZADO (A)")
println("-" * 70)
for (i, p) in enumerate(polos)
    println("  Polo $(lpad(i,2)):  $(round(real(p), digits=4))  +  $(round(imag(p), digits=4))im")
end

instaveis = count(real(polos) .> 1e-6)
println("-" * 70)
if instaveis > 0
    println("  ⚠️  Sistema possui $instaveis polo(s) instável(is).")
else
    println("  ✅ Sistema estável.")
end
println("=" * 70)