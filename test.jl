using OrdinaryDiffEq
using PlotlyJS 

include("Controle/Funções Auxiliares/CoversãoPoly2EDO.jl")
include("Controle/Funções Auxiliares/ResoluçãoEDO.jl")
include("Controle/Funções Auxiliares/Visualização.jl")

# 1. Suas variáveis e sistema
@polyvar x1 x2
vars = [x1, x2]

f = [2*x1^3 + x1^2*x2 - 6*x1*x2^2 + 5*x2^3, 0.0]
g = [0.0, 1.0]
u = 1.0582933631271871e-29 + 0.0012313160167337946*x2 - 0.0018180597277442838*x1 + 2.0678052773368813e-12*x2^2 - 3.722042803877773e-12*x1*x2 + 2.2932506152804264e-12*x1^2 - 3.558277299249908*x2^3 + 1.5723270323261815*x1*x2^2 + 3.7987616032273896*x1^2*x2 - 3.957236453645706*x1^3
u1 = -3.6345*x1^3 + 4.4439*x1^2*x2 - 7.5113*x1*x2^2 - 3.5452*x2^3;

meu_sistema_numerico1! = extrairSistemaEDO(f, vars)
meu_sistema_numerico2! = extrairSistemaEDO(f, vars, g_x=g, u_x=u)
meu_sistema_numerico3! = extrairSistemaEDO(f, vars, g_x=g, u_x=u1)

x0    = [0.5, -0.5]
tspan = (0.0, 1000.0)

solucao = resolverSistema(meu_sistema_numerico1!, x0, tspan)
solucao1 = resolverSistema(meu_sistema_numerico2!, x0, tspan)
solucao2 = resolverSistema(meu_sistema_numerico3!, x0, tspan)

plotarNoTempo([solucao, solucao1, solucao2])
# plotarNoTempo(solucao1)

plotarRetratoFase([solucao, solucao1, solucao2])
# plotarRetratoFase(solucao1)