using OrdinaryDiffEq
using Plots

# 1. Definir a função do sistema (du/dt = f(u, p, t))
# u[1] é a posição (x), u[2] é a velocidade (v)
# p são os parâmetros [k, m]
function massa_mola!(dx, x, constantes, t)
    k = constantes.k; b = constantes.b; m = constantes.m; 

    dx[1] = x[2]           # dx/dt = v
    dx[2] = -(b/m) * x[2] - (k/m) * x[1]  # dv/dt = -(k/m)x
end

# 2. Condições Iniciais e Parâmetros
u0 = [1.0, 0.0]      
tspan = (0.0, 20.0)   
constantes = (k=10.0, b=1, m=1)    

# 3. Definir e Resolver o Problema
prob = ODEProblem(massa_mola!, u0, tspan, constantes)
sol  = solve(prob)

# 4. Plotar os resultados
plot(sol, idxs=(0, 1), title="Sistema Massa-Mola",
     xlabel="Tempo (t)", ylabel="Posição (x)", label="Posição")