using JuMP
using SumOfSquares
using DynamicPolynomials
using COSMO
using OrdinaryDiffEq
using LinearAlgebra

function encontrarControleSOS(f, g, V, vars)
    # model = SOSModel(Mosek.Optimizer)
    model = SOSModel(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "eps_abs" => 1e-5))
    set_silent(model) 

    # Define os monômios
    Z_u = monomials(vars, 0:3)
    Z_λ = monomials(vars, 2:6)

    # Cria as variáveis polinomiais de decisão
    @variable(model, u, Poly(Z_u))
    @variable(model, λ, Poly(Z_λ)) 

    # Calcula o gradiente e a restrição principal
    grad_V = differentiate(V, vars)
    expr   = -dot(grad_V, f + g*u) + λ

    # Restrição de positividade SOS
    @constraint(model, expr in SOSCone())

    # Resolve o problema
    optimize!(model)
    status = termination_status(model)

    # No JuMP, MOI.OPTIMAL garante que o problema foi resolvido com sucesso (substitui a checagem de feasratio/pinf)
    if status == MOI.OPTIMAL
        # println("Solução encontrada com sucesso!")
        return value(u), value(λ), true
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

    # Restrição 1: V(x) - 1e-2*(x1^6 + x2^6) >= 0
    @constraint(model, V - 1e-2*(vars[1]^6 + vars[2]^6) in SOSCone())

    # Restrição 2: -dV/dx * f_cl >= 0
    grad_V = differentiate(V, vars)
    expr   = -dot(grad_V, f + g*u)
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