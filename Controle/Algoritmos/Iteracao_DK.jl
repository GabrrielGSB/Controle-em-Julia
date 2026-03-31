using JuMP
using SumOfSquares
using DynamicPolynomials
# using MosekTools 
using COSMO
using DifferentialEquations
using Plots
using LinearAlgebra
plotlyjs()

# =========================================================================
# Funções Auxiliares (SOS)
# =========================================================================

function encontrarControleSOS(f, g, V, vars)
    # model = SOSModel(Mosek.Optimizer)
    model = SOSModel(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "eps_abs" => 1e-5))
    set_silent(model) 

    # Define os monômios
    Z_u = monomials(vars, 0:3)
    Z_lamb = monomials(vars, 2:6)

    # Cria as variáveis polinomiais de decisão
    @variable(model, u, Poly(Z_u))
    @variable(model, lamb, Poly(Z_lamb)) 

    # Calcula o gradiente e a restrição principal
    grad_V = differentiate(V, vars)
    expr = -dot(grad_V, f + g*u) + lamb

    # Restrição de positividade SOS
    @constraint(model, expr in SOSCone())

    # Resolve o problema
    optimize!(model)
    status = termination_status(model)

    # No JuMP, MOI.OPTIMAL garante que o problema foi resolvido com sucesso (substitui a checagem de feasratio/pinf)
    if status == MOI.OPTIMAL
        # println("Solução encontrada com sucesso!")
        return value(u), value(lamb), true
    else
        @warn "O solver falhou ao encontrar u(x) (Status: $status)"
        return nothing, nothing, false
    end
end

function encontrarLyapunovSOS(f, g, u, vars)
    # model = SOSModel(Mosek.Optimizer)
    model = SOSModel(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "eps_abs" => 1e-5))
    set_silent(model)

    Z_V = monomials(vars, 2:6)
    @variable(model, V, Poly(Z_V))

    f_cl = f + g*u

    # Restrição 1: V(x) - 1e-2*(x1^6 + x2^6) >= 0
    @constraint(model, V - 1e-2*(vars[1]^6 + vars[2]^6) in SOSCone())

    # Restrição 2: -dV/dx * f_cl >= 0
    grad_V = differentiate(V, vars)
    expr = -dot(grad_V, f_cl)
    @constraint(model, expr in SOSCone())

    optimize!(model)
    status = termination_status(model)

    if status == MOI.OPTIMAL
        # println("Função de Lyapunov V(x) encontrada com sucesso!")
        return value(V), true
    else
        @warn "O solver falhou ao encontrar V(x) (Status: $status)"
        return nothing, false
    end
end

# =========================================================================
# Função de Simulação e Plotagem
# =========================================================================

function plotar_comparacao(f, g, u, vars, x0, t_final)
    f_fechada = f + g*u

    # Converte os polinômios simbólicos para funções numéricas eficientes
    sys_aberta!(dx, x, p, t) = begin
        dx[1] = f[1](vars[1] => x[1], vars[2] => x[2])
        dx[2] = f[2](vars[1] => x[1], vars[2] => x[2])
    end

    sys_fechada!(dx, x, p, t) = begin
        dx[1] = f_fechada[1](vars[1] => x[1], vars[2] => x[2])
        dx[2] = f_fechada[2](vars[1] => x[1], vars[2] => x[2])
    end

    # Callback de segurança para evitar que a simulação vá para o infinito (norma > 100)
    condicao_explosao(x, t, integrator) = norm(x) > 100
    trava_explosao = DiscreteCallback(condicao_explosao, terminate!)

    tspan = (0.0, t_final)
    
    println("\nSimulando malha aberta...")
    prob_ab = ODEProblem(sys_aberta!, x0, tspan)
    sol_ab  = solve(prob_ab, Tsit5(), callback=trava_explosao)

    println("Simulando malha fechada...")
    prob_fe = ODEProblem(sys_fechada!, x0, tspan)
    sol_fe  = solve(prob_fe, Tsit5(), callback=trava_explosao)

    # Extraindo dados para os plots
    t_ab, x1_ab, x2_ab = sol_ab.t, [u[1] for u in sol_ab.u], [u[2] for u in sol_ab.u]
    t_fe, x1_fe, x2_fe = sol_fe.t, [u[1] for u in sol_fe.u], [u[2] for u in sol_fe.u]

    # Layout e Plotagem usando Plots.jl
    l = @layout [a c; b c]
    
    p1 = plot(t_ab, x1_ab, label="Malha Aberta", line=(:red, :dash), lw=2)
    plot!(p1, t_fe, x1_fe, label="Malha Fechada", line=(:blue, :solid), lw=2, title="Evolução de x_1", xlabel="Tempo (s)", ylabel="x_1")

    p2 = plot(t_ab, x2_ab, label="Malha Aberta", line=(:red, :dash), lw=2)
    plot!(p2, t_fe, x2_fe, label="Malha Fechada", line=(:blue, :solid), lw=2, title="Evolução de x_2", xlabel="Tempo (s)", ylabel="x_2")

    p3 = plot(x1_ab, x2_ab, label="Trajetória M. Aberta", line=(:red, :dash), lw=2)
    plot!(p3, x1_fe, x2_fe, label="Trajetória M. Fechada", line=(:blue, :solid), lw=2)
    scatter!(p3, [x0[1]], [x0[2]], color=:green, markersize=6, label="Início")
    scatter!(p3, [0.0], [0.0], color=:black, shape=:xcross, markersize=8, label="Origem")
    plot!(p3, title="Plano de Fase (x_1 vs x_2)", xlabel="x_1", ylabel="x_2", aspect_ratio=:equal)

    display(plot(p1, p2, p3, layout=l, size=(900, 600)))
end

# =========================================================================
# Main: Algoritmo de Iteração D-K
# =========================================================================

# 1. Definição das Variáveis e do Sistema
@polyvar x1 x2
vars = [x1, x2]

# Dinâmica do sistema (f) e vetor de entrada (g)
f = [2*x1^3 + x1^2*x2 - 6*x1*x2^2 + 5*x2^3, 0.0]
g = [0.0, 1.0]

# 2. Configurações da Iteração
max_iter = 1 
iter = 1

# 3. Estimativa V inicial (Chute pelo método do controle virtual)
V_atual = 2*x1^2 + 2*x1*x2 + x2^2 

println("=========================================")
println("Iniciando Iteração D-K")
println("Chute Inicial V0: ", V_atual)
println("=========================================\n")

u_final = nothing
V_final = V_atual

# 4. Loop Principal
while iter <= max_iter
    println("--- Iteração $iter ---")
   
    # Passo 1
    println("  > Buscando controle u(x)...")
    u_atual, lamb_atual, success_u = encontrarControleSOS(f, g, V_atual, vars)
    
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
    println("  u = ", u_atual, "\n")
    
    global iter += 1
end

# =========================================================================
# Resumo Final
# =========================================================================
println("\n=========================================")
println("Fim do Algoritmo D-K (Iterações: $(iter - 1))")
println("=========================================")

if u_final === nothing
    println("O algoritmo não conseguiu estabilizar o sistema já na 1ª iteração.")
    println("Tente um chute inicial diferente para V0(x).")
else
    println("Melhor controle u(x) viável encontrado:")
    println(u_final, "\n")
    
    println("Melhor função V(x) viável encontrada:")
    println(V_final, "\n")
    
    println("Gerando gráficos de comparação...")
    condicao_inicial = [1.0, -1.0]
    tempo_simulacao = 200.0 
    
    plotar_comparacao(f, g, u_final, vars, condicao_inicial, tempo_simulacao)
end