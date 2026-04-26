# Controle/Ferramentas/ResoluçãoEDO.jl

using OrdinaryDiffEq
using DiffEqCallbacks

function resolverSistema(sistema, x0, intervaloTempo, p=sistema; 
                         resolucao=0.01,
                         salvar_controle=false)   # ← nova opção

    x0 = x0 isa Vector{Float64} ? x0 : collect(Float64, x0)

    # Histórico de controle — preenchido pelo callback
    historico_u = Vector{Vector{Float64}}()
    historico_t = Vector{Float64}()

    # Callback: chamado pelo solver a cada passo salvo
    cb = if salvar_controle 
            SavingCallback((x, t, integrador) -> begin
                    n = integrador.p.dim_planta
                    x_planta = x[1:n]
                    x_ctrl   = x[n+1:end]
                    u = calcularSaida(integrador.p.controlador,
                                    x_planta, x_ctrl,
                                    integrador.p.referencia, t)
                    # Normaliza para sempre ser vetor
                    u_vec = u isa AbstractVector ? u : [u]
                    push!(historico_t, t)
                    push!(historico_u, copy(u_vec))
                    nothing
                end,
                SavedValues(Float64, Nothing)
            )
         else
             nothing
         end

    problema = if p === nothing
        ODEProblem(sistema, x0, intervaloTempo)
    else
        ODEProblem(sistema, x0, intervaloTempo, p)
    end

    solucao = if cb === nothing
        solve(problema, saveat=resolucao)
    else
        solve(problema, saveat=resolucao, callback=cb)
    end

    # Retorna tupla nomeada quando histórico foi pedido
    if salvar_controle
        return (solucao=solucao, 
                t_u=historico_t, 
                u=historico_u)
    else
        return solucao
    end
end