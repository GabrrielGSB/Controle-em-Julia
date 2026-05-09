#=
SISTEMA DINÂMICO: PÊNDULO INVERTIDO NO CARRINHO (CART-POLE)
-----------------------------------------------------------------
DESCRIÇÃO: Representa um carrinho de massa 'M' que se move no eixo X 
sob a ação de uma força 'u', equilibrando uma haste de comprimento 'l' 
e massa 'm' na ponta.

ESTADOS:
x[1] = Posição do carrinho (m)
x[2] = Velocidade do carrinho (m/s)
x[3] = Ângulo do pêndulo (rad) -> 0 é para CIMA (invertido)
x[4] = Velocidade angular do pêndulo (rad/s)
=#

using ProgressMeter
import Plots
include("../Abstrações/Interfaces.jl")

#=========================================================================
STRUCT DA PLANTA
=========================================================================#

struct PenduloInvertido <: Planta
    M               ::Float64  # Massa do carrinho (kg)
    m               ::Float64  # Massa do pêndulo (kg)
    l               ::Float64  # Comprimento da haste (m)
    g               ::Float64  # Aceleração da gravidade (m/s²)
    b               ::Float64  # Atrito viscoso do carrinho (N.s/m)
    estadosIniciais ::Vector{Float64}
    numEstados      ::Int
    variaveisEstado ::Tuple{String, String, String, String}
    nomeSistema     ::String
    dinamica!       ::Function # NOVIDADE: Adicionando o campo para armazenar a função
end

"""
Construtor amigável com valores padrão clássicos de benchmark.
"""
function PenduloInvertido(; 
        M = 1.0, 
        m = 0.1, 
        l = 0.5, 
        g = 9.81, 
        b = 0.1, 
        estadosIniciais = [0.0, 0.0, 0.1, 0.0] # Começa com pequeno desvio angular
    )
    
    return PenduloInvertido(
        M, m, l, g, b, 
        estadosIniciais, 
        4, 
        ("Posição Carro (m)", "Velocidade Carro (m/s)", "Ângulo (rad)", "Velocidade Angular (rad/s)"), 
        "Pêndulo Invertido no Carrinho",
        pendulo_invertido! # NOVIDADE: Vinculando a função à instância da struct
    )
end

#=========================================================================
DINÂMICA DA PLANTA 
=========================================================================#

function pendulo_invertido!(dx, x, p, t; u=0.0)
    # Extração de parâmetros para limpeza das equações
    M, m, l, g, b = p.M, p.m, p.l, p.g, p.b
    
    # Extração dos estados
    x_c = x[1]
    v_c = x[2]
    θ   = x[3]
    ω   = x[4]
    
    # Termos trigonométricos recorrentes
    senθ = sin(θ)
    cosθ = cos(θ)
    
    # Denominador comum (matriz de inércia)
    # O sistema é não-linear, a inércia efetiva muda com o ângulo
    D = M + m * (senθ^2)
    
    # Dinâmica translacional do carrinho
    dx[1] = v_c
    dx[2] = (u - b*v_c + m*l*(ω^2)*senθ - m*g*senθ*cosθ) / D
    
    # Dinâmica rotacional do pêndulo
    dx[3] = ω
    dx[4] = ((M + m)*g*senθ - cosθ*(u - b*v_c + m*l*(ω^2)*senθ)) / (l * D)
end



#=========================================================================
ANIMAÇÃO (Opcional - Estrutura Base)
=========================================================================#

function gerarAnimacao(solucao, p::PenduloInvertido, fps)
    # A implementação da animação gráfica exigirá renderizar um retângulo 
    # (carrinho) e uma linha com um círculo na ponta (pêndulo).
    println("Animação do Cart-Pole em construção...")
end