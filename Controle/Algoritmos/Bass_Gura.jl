
include("../Ferramentas/Polinômios.jl")
include("../Ferramentas/LinearizarSistema.jl")
include("../Ferramentas/PropriedadesSistema.jl")

"""
    bass_gura(A, B, polos_desejados)

    Aplica a Fórmula de Bass-Gura para encontrar o vetor de ganhos K 
    que aloca os polos de um sistema linear de única entrada (SISO).
"""
function bass_gura(A::Matrix{Float64}, B::Vector{Float64}, polos_desejados::Vector{<:Number})
    n = size(A, 1)
    
    if (!controlavel(A,B)) error("O sistema não é controlável!")
    else                   P = matrizControlabilidade(A,B) end
    
    lambdas_A = eigvals(A)
    a_coefs  = coeficientesPolinomiais(lambdas_A)
    
    alpha_coefs = coeficientesPolinomiais(polos_desejados)
    
    # d) Matriz W (Toeplitz superior modificada com os coeficientes de malha aberta)
    W = zeros(n, n)
    for i in 1:n
        for j in 1:n
            k = i + j
            if k <= n + 1
                W[i, j] = a_coefs[k] 
            else
                W[i, j] = 0.0
            end
        end
    end
    
    diff_coefs = alpha_coefs[1:n] - a_coefs[1:n]
    
    # Aplicação da Fórmula de Bass-Gura: K = Δa^T * (P * W)^-1
    K_transposto = diff_coefs' * inv(P * W)
    
    return K_transposto' 
end
