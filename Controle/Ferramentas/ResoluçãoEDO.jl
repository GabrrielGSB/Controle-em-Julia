using OrdinaryDiffEq
using DiffEqCallbacks

# =========================================================================
# SOLVER
    function resolverSistema(sistema, x0, intervaloTempo, p=sistema; 
                         resolucao=0.01,
                         salvar_controle=false)

    x0 = x0 isa Vector{Float64} ? x0 : collect(Float64, x0)

    problema = ODEProblem(sistema, x0, intervaloTempo, p)
    solucao  = solve(problema, saveat=resolucao)

    if salvar_controle
        n         = p.dim_planta
        historico_u = Vector{Vector{Float64}}()
        historico_t = solucao.t   

        for (i, t) in enumerate(solucao.t)
            x       = solucao.u[i]
            x_planta = x[1:n]
            x_ctrl   = x[n+1:end]
            u = calcularSaida(p.controlador, x_planta, x_ctrl, p.referencia, t)
            u_vec = u isa AbstractVector ? u : [u]
            push!(historico_u, copy(u_vec))
        end

        return (solucao=solucao, t_u=historico_t, u=historico_u)
    else
        return solucao
    end
end
# =========================================================================