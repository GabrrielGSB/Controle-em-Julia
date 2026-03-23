using OrdinaryDiffEq
using PlotlyJS

function simularSistemaMalhaAberta(; sistema, estadosIniciais, parametros, tempoSimulacao, 
                                     visualizacao=true, eixosVisu=(0, 1), animar = false, animarFPS = 30, CSV=false)

    estadosIniciais = collect(values(estadosIniciais))

    problema = ODEProblem(sistema, estadosIniciais, tempoSimulacao, parametros)
    solucao  = solve(problema)

    if visualizacao
        plotlyjs() 
        display(Plots.plot(solucao, 
                   idxs   = eixosVisu, 
                   title  = parametros.nomeSistema,
                   xlabel =  eixosVisu[1] == 0 ? "Tempo (t)" : parametros.variaveisEstado[eixosVisu[1]], 
                   ylabel =  eixosVisu[2] == 0 ? "Tempo (t)" : parametros.variaveisEstado[eixosVisu[2]], 
                   label  = "Comportamento do sistema",
                   lw     = 2,       # Largura da linha
                   marker = :circle, # Adiciona pontos para você passar o mouse por cima
                   ms     = 2))       # Tamanho do marcador
    end 

    if animar   
        animacaoPendulo(solucao = solucao, 
                        params  = parametros, 
                        fps     = animarFPS)
    end
end