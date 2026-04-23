# Controle/Ferramentas/ResoluçãoEDO.jl

using OrdinaryDiffEq

"""
    resolverSistema(sistema, x0, intervaloTempo, p; resolucao) -> solucao

    Resolve uma EDO genérica e retorna o objeto de solução.

    Argumentos:
        - sistema        : função f(dx,x,p,t) OU struct callable (ex: MalhaFechada)
        - x0             : vetor de condições iniciais
        - intervaloTempo : tupla (t_inicial, t_final)
        - p              : parâmetros externos — ignorado se sistema for callable struct
        - resolucao      : intervalo de salvamento dos pontos (padrão 0.01s)
"""
function resolverSistema(sistema, x0, intervaloTempo, p=nothing; resolucao=0.01)
    x0 = x0 isa Vector{Float64} ? x0 : collect(Float64, x0)

    problema = if p === nothing
        ODEProblem(sistema, x0, intervaloTempo)
    else
        ODEProblem(sistema, x0, intervaloTempo, p)
    end

    solucao = solve(problema, saveat=resolucao)

    return solucao
end

using OrdinaryDiffEq

function resolverSistema(sistema, x0, intervaloTempo, p=nothing; 
                         resolucao=0.01,
                         salvar_controle=false)   # ← nova opção

    x0 = x0 isa Vector{Float64} ? x0 : collect(Float64, x0)

    # Histórico de controle — preenchido pelo callback
    historico_u = Vector{Vector{Float64}}()
    historico_t = Vector{Float64}()

    # Callback: chamado pelo solver a cada passo salvo
    cb = if salvar_controle && p isa MalhaFechada
        SavingCallback(
            (x, t, integrador) -> begin
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