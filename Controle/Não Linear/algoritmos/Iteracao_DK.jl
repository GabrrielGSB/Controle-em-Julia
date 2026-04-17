include("../../Funções Auxiliares/SOS.jl")

# 1. Definição das Variáveis e do Sistema
@polyvar x1 x2
vars = [x1, x2]

# Dinâmica do sistema (f) e vetor de entrada (g)
f = [2*x1^3 + x1^2*x2 - 6*x1*x2^2 + 5*x2^3, 0.0]
g = [0.0, 1.0]

# 2. Configurações da Iteração
max_iter = 1 

# 3. Estimativa V inicial 
V_atual = 2*x1^2 + 2*x1*x2 + x2^2 

println("=========================================")
println("Iniciando Iteração D-K")
println("Chute Inicial V0: ", V_atual)
println("=========================================\n")

u_final = nothing
V_final = V_atual

# 4. Loop Principal
iter = 1
while iter <= max_iter
    println("--- Iteração $iter ---")
   
    # Passo 1
    println("  > Buscando controle u(x)...")
    u_atual, λ_atual, success_u = encontrarControleSOS(f, g, V_atual, vars)
    
    if !success_u
        println("  [!] Solver falhou ao achar u(x). Parando iterações.")
        break
    end
    
    # Passo 2
    println("  > Buscando nova função V(x)...")
    V_novo, success_V = encontrarLyapunovSOS(f, g, u_atual, vars)
    
    if !success_V
        println("  [!] Solver falhou ao achar V(x). Parando iterações.")
        break
    end
   
    # Passo 3
    global V_atual = V_novo 
    global V_final = V_novo
    global u_final = u_atual
    
    println("  [OK] Novo par (u, V) viável encontrado!")
    # O Julia já formata o polinômio maravilhosamente bem no console
    # println("  u = ", u_atual, "\n")
    
    global iter += 1
end

println("\n=========================================")
println("Fim do Algoritmo D-K (Iterações: $(iter - 1))")
println("=========================================")

if u_final === nothing
    println("O algoritmo não conseguiu estabilizar o sistema já na 1ª iteração.")
    println("Tente um chute inicial diferente para V0(x).")
else
    println("Melhor controle u(x) viável encontrado:")
    println(u_final, "\n")
    
    # println("Melhor função V(x) viável encontrada:")
    # println(V_final, "\n")
end