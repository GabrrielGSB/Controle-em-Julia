using OrdinaryDiffEq
using Plots

include("Controle/SistemasLTI.jl")

# Função do Pêndulo Simples
function pendulo!(dx, x, p, t)
    g, L, b, m = p
    dx[1] = x[2]
    dx[2] = -(g/L)*sin(x[1]) - (b/(m*L^2))*x[2]
end

# Parâmetros e Solução
p = [9.81, 1.0, 0.1, 1.0] # g, L, b, m
u0 = [pi/4, 0.0]           # 45 graus inicial
tspan = (0.0, 20.0)
prob = ODEProblem(pendulo!, u0, tspan, p)
sol = solve(prob) # saveat define a "fluidez" da animação (FPS)

animacaoPendulo(solucao=sol, params=p, fps=30)