using OrdinaryDiffEq

"""
    Resolve uma EDO genérica e retorna o objeto de solução.
    - sistema: A função que define dx/dt (a dinâmica).
    - x0: Vetor de condições iniciais.
    - tspan: Intervalo de tempo
    - p: Parâmetros do sistema (opcional).
    - resolucao: Forçar o salvamento de pontos a cada 't' tempo (padrão->0.01s).
"""
function resolverSistema(sistema, condicoesIniciais, intervaloTempo, p=nothing; 
                         resolucao=0.01)

    condicoesIniciais = collect(values(condicoesIniciais))

    problema = ODEProblem(sistema, condicoesIniciais, intervaloTempo, p)
    solucao  = solve(problema, saveat=resolucao)

    return solucao
end