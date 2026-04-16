using OrdinaryDiffEq
using PlotlyJS

function simularSistemaMalhaAberta(; sistema, estadosIniciais, parametros, tempoSimulacao, 
                                     visualizacao=true, eixosVisu=(0, 1), 
                                     visualizarTodosEstados = false,
                                     animar = false, animarFPS = 30, CSV=false)

    problema = ODEProblem(sistema, estadosIniciais, tempoSimulacao, parametros)
    solucao  = solve(problema)

    if visualizarTodosEstados
        plotlyjs() 
        
        # Converte a tupla de nomes em uma matriz 1xN para aplicar os nomes corretos na legenda
        num_estados = length(parametros.variaveisEstado)
        nomes_matriz = reshape(collect(parametros.variaveisEstado), 1, num_estados)

        display(Plots.plot(solucao, 
                   title  = "$(parametros.nomeSistema) - Todos os Estados",
                   xlabel = "Tempo (t)", 
                   ylabel = "Amplitude", 
                   label  = nomes_matriz,
                   lw     = 4,       # Largura da linha
                   marker = :circle, # Adiciona pontos interativos
                   ms     = 0.5))      # Tamanho do marcador
    end

    if visualizacao
        plotlyjs() 
        display(Plots.plot(solucao, 
                   idxs   = eixosVisu, 
                   title  = parametros.nomeSistema,
                   xlabel =  eixosVisu[1] == 0 ? "Tempo (t)" : parametros.variaveisEstado[eixosVisu[1]], 
                   ylabel =  eixosVisu[2] == 0 ? "Tempo (t)" : parametros.variaveisEstado[eixosVisu[2]], 
                   label  = "Comportamento",
                   lw     = 2,       # Largura da linha
                   marker = :circle, # Adiciona pontos para você passar o mouse por cima
                   ms     = 2))       # Tamanho do marcador
    end 

    if animar   
        gerarAnimacao(solucao, parametros, animarFPS)
    end
end