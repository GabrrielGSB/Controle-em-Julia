using OrdinaryDiffEq
using PlotlyJS
using Plots

# Ativa o backend interativo
plotlyjs() 

function massa_mola!(dx, x, constantes, t)
    k, b, m = constantes.k, constantes.b, constantes.m
    dx[1] = x[2]
    dx[2] = -(b/m) * x[2] - (k/m) * x[1]
end

u0 = [1.0, 0.0]      
tspan = (0.0, 20.0)   
constantes = (k=10.0, b=1.0, m=1.0)    

prob = ODEProblem(massa_mola!, u0, tspan, constantes)
sol  = solve(prob)

# O comando de plot permanece quase igual, 
# mas agora ele terá ferramentas de zoom e hover automáticas
Plots.plot(sol, 
           idxs=(0, 1), 
           title="Sistema Massa-Mola Interativo",
           xlabel="Tempo (t)", 
           ylabel="Posição (x)", 
           label="Posição (m)",
           lw=2,              # Largura da linha
           marker=:circle,    # Adiciona pontos para você passar o mouse por cima
           ms=2,              # Tamanho do marcador
          )