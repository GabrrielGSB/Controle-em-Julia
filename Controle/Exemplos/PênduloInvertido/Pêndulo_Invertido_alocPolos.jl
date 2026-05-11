using ForwardDiff
using LinearAlgebra

include("../../Ferramentas/LinearizarSistema.jl")
include("../../SistemasLTI/PênduloInvertido.jl")


# =======================================================================
# 2. IMPLEMENTAÇÃO DA FÓRMULA DE BASS-GURA
# =======================================================================
"""
    bass_gura(A, B, polos_desejados)

    Aplica a Fórmula de Bass-Gura para encontrar o vetor de ganhos K 
    que aloca os polos de um sistema linear de única entrada (SISO).
"""
function bass_gura(A::Matrix{Float64}, B::Vector{Float64}, polos_desejados::Vector{<:Number})
    n = size(A, 1)
    
    # a) Matriz de Controlabilidade P = [B, AB, A²B, ..., A^(n-1)B]
    P = zeros(n, n)
    P[:, 1] = B
    for i in 2:n
        P[:, i] = A * P[:, i-1]
    end
    
    if abs(det(P)) < 1e-8
        error("O sistema não é controlável! Determinante da Matriz de Controlabilidade é zero.")
    end
    
    # Função interna: Calcula os coeficientes do polinômio a partir das raízes
    # Retorna [a_0, a_1, ..., a_n] onde a_n = 1
    function poly_coefs(raizes)
        p = zeros(ComplexF64, n + 1)
        p[1] = 1.0 # Polinômio inicial: P(s) = 1 (s^0)
        p_temp = zeros(ComplexF64, n + 1)
        
        for r in raizes
            p_temp .= 0.0
            # Multiplicar por s (desloca os coeficientes)
            for i in 1:n
                p_temp[i+1] += p[i]
            end
            # Multiplicar por -r
            for i in 1:n
                p_temp[i] -= r * p[i]
            end
            p .= p_temp
        end
        return real.(p) # Como sistemas físicos têm raízes conjugadas, a parte imaginária zera
    end
    
    # b) Coeficientes de Malha Aberta (a_i)
    lambdas_A = eigvals(A)
    a_coefs = poly_coefs(lambdas_A)
    
    # c) Coeficientes de Malha Fechada Desejados (alpha_i)
    alpha_coefs = poly_coefs(polos_desejados)
    
    # d) Matriz W (Toeplitz superior modificada com os coeficientes de malha aberta)
    W = zeros(n, n)
    for i in 1:n
        for j in 1:n
            k = i + j
            if k <= n + 1
                W[i, j] = a_coefs[k] # Lembre-se: índice k = 1 corresponde a a_0, etc.
            else
                W[i, j] = 0.0
            end
        end
    end
    
    # e) Diferença dos coeficientes: [alpha_0 - a_0, ..., alpha_{n-1} - a_{n-1}]
    diff_coefs = alpha_coefs[1:n] - a_coefs[1:n]
    
    # f) Aplicação da Fórmula de Bass-Gura: K = Δa^T * (P * W)^-1
    K_transposto = diff_coefs' * inv(P * W)
    
    return K_transposto' # Retornamos um vetor (n x 1)
end

# =======================================================================
# 3. DEFINIÇÃO DA PLANTA (PÊNDULO INVERTIDO)
# =======================================================================
# Adaptando a struct do seu frame work para simulação local
struct PenduloInvertido
    M::Float64   # Massa do carrinho
    m::Float64   # Massa do pêndulo
    l::Float64   # Comprimento da haste
    g::Float64   # Gravidade
    b::Float64   # Atrito
    numEstados::Int
    dinamica!::Function
end

function pendulo_dinamica!(dx, x, p, t; u=0.0)
    # Estados: x1 = p, x2 = p_ponto, x3 = theta (0 para cima), x4 = theta_ponto
    M, m, l, g, b = p.M, p.m, p.l, p.g, p.b
    
    x1, x2, x3, x4 = x[1], x[2], x[3], x[4]
    
    # Dinâmica Não-Linear
    sen_t = sin(x3)
    cos_t = cos(x3)
    denominador = M + m - m * cos_t^2
    
    dx[1] = x2
    dx[2] = (u + m * l * x4^2 * sen_t - m * g * sen_t * cos_t - b * x2) / denominador
    dx[3] = x4
    dx[4] = ((M + m) * g * sen_t - cos_t * (u + m * l * x4^2 * sen_t - b * x2)) / (l * denominador)
end

# =======================================================================
# 4. EXECUÇÃO DO PROJETO
# =======================================================================
# Instancia a planta
planta_pendulo = PenduloInvertido(1.0, 0.1, 0.5, 9.81, 0.1, 4, pendulo_dinamica!)

# Define o ponto de operação (Pêndulo equilibrado para cima)
x_eq = [0.0, 0.0, 0.0, 0.0]
u_eq = [0.0]

# Extrai o modelo linear
println("1. Linearizando o Sistema...")
A, B_mat, C, D = linearizar_sistema(planta_pendulo, x_eq, u_eq)
B = B_mat[:, 1] # Converte para vetor coluna unidimensional exigido por Bass-Gura

display(A)
display(B)

# Define onde queremos os polos do sistema em malha fechada
# (Exemplo: Polos reais e conjugados estáveis)
polos_alvo = [-2.0 + 1.0im, -2.0 - 1.0im, -3.0, -4.0]

println("\n2. Calculando Ganhos via Bass-Gura...")
K = bass_gura(A, B, polos_alvo)

println("\nVetor de Ganhos K encontrado:")
display(K')

println("\n3. Validando Polos de Malha Fechada (A - B*K):")
A_fechada = A - B * K'
polos_obtidos = eigvals(A_fechada)
display(polos_obtidos)