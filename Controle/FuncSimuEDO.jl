using OrdinaryDiffEq
using PlotlyJS
using Plots

function simularSistemaMalhaAberta(; sistema, estadosIniciais, parametros, tempoSimulacao, visualizacao=true, CSV=false)
    problema = ODEProblem(sistema, estadosIniciais, tempoSimulacao, parametros)
    solucao = solve(problema)

    if visualizacao
        Plots.plot(solucao, 
                   idxs   = (0, 1), 
                   title  = "Sistema",
                   xlabel = "Tempo (t)", 
                   ylabel = "Posição (rad)", 
                   label  = "Posição (rad)",
                   lw     = 2,       # Largura da linha
                   marker = :circle, # Adiciona pontos para você passar o mouse por cima
                   ms     = 2)       # Tamanho do marcador     
    end 
end