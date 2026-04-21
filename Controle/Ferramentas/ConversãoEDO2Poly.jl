using DynamicPolynomials

import Base: sin, cos, tan

sin(p::AbstractPolynomialLike) = p - (p^3)/6 
cos(p::AbstractPolynomialLike) = 1 - (p^2)/2 
tan(p::AbstractPolynomialLike) = p + (p^3)/3


"""
    Recebe uma função de EDO numérica padrão e retorna os polinômios f(x), g(x) 
    e as variáveis simbólicas compatíveis com o controle SOS.
"""
function extrairSistemaAfim(sistema!, params, dim_x; g=[0.0, 1.0])
    # 1. Cria as variáveis simbólicas usando DynamicPolynomials
    @polyvar x[1:dim_x]
    
    # 2. Prepara um vetor vazio flexível para receber as equações
    f = Vector{Any}(undef, dim_x)
    
    # 3. Roda a sua função numérica! 
    sistema!(f, x, params, 0.0)
    
    # 4. Define a matriz de entrada do controle g(x).
    g = g
    
    return f, g, x
end