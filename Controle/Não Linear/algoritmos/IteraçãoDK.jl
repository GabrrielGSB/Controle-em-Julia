include("../../Funções Auxiliares/SOS.jl")

function iteracaoDK(f, g, vars, V_inicial; max_iter=4, nome_subsistema="GENÉRICO")
    println("===================================================================")
    println("Iniciando D-K Iteration - Subsistema de $nome_subsistema")
    println("===================================================================")

    # Inicialização de variáveis locais
    V_atual = V_inicial
    u_final = nothing
    V_final = V_atual
    iter = 1

    while iter <= max_iter
        println("--- Iteração $iter ---")
       
        # Passo 1: Otimizar o Controle (u)
        u_atual, lamb_atual, success_u = encontrarControleSOS(f, g, V_atual, vars)
        
        if !success_u
            println("  [!] Solver falhou ao achar u(x). Parando iterações.")
            break
        end
        
        # Passo 2: Otimizar a Função de Lyapunov (V)
        V_novo, success_V = encontrarLyapunovSOS(f, g, u_atual, vars)
        
        if !success_V
            println("  [!] Solver falhou ao achar V(x). Parando iterações.")
            break
        end
       
        # Atualização das variáveis no escopo da função (sem precisar de 'global')
        V_atual = V_novo 
        V_final = V_novo
        u_final = u_atual
        
        println("--> Controlador u(x) estabiliza o sistema!")
        
        println("  u(x) = ", round.(u_atual, digits=4), "\n") 
        
        iter += 1
    end

    println("===================================================================\n\n")
    if u_final === nothing
        println("O algoritmo falhou na estabilização de $nome_subsistema.")
    else
        # println("SUCESSO! O controlador de $nome_subsistema foi estabilizado.")
        # println("\nEquação de Controle Virtual encontrada:")
        # println("U_barra = ", round.(u_final, digits=4))
    end
    
    # return u_final
end