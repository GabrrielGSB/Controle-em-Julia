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