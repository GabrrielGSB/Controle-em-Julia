# =======================================================================
# ARQUIVO: Controle/Ferramentas/Polinomios.jl
# =======================================================================

"""
    poly_coefs(raizes::Vector{<:Number}) -> Vector{Float64}

Calcula iterativamente os coeficientes reais de um polinômio a partir de 
um conjunto de raízes (reais ou complexas).

# Matemática
Dado um vetor de raízes `r = [r₁, r₂, ..., rₙ]`, esta função computa os 
coeficientes `a` do polinômio característico expandido:
P(s) = (s - r₁)(s - r₂)...(s - rₙ) = a₀ + a₁s + a₂s² + ... + aₙsⁿ

A expansão é feita multiplicando o polinômio atual por (s - rᵢ) a cada iteração, 
evitando a necessidade de pacotes externos simbólicos.

# Argumentos
- `raizes`: Um vetor contendo as raízes desejadas para o polinômio. Pode conter 
números Reais ou Complexos. Sistemas físicos devem usar raízes complexas em 
pares conjugados.

# Retorno
- `Vector{Float64}`: Um vetor de tamanho `n+1` contendo os coeficientes ordenados 
do menor grau (termo independente `a₀`) para o maior grau (`aₙ = 1.0`). A 
parte imaginária é descartada sob a premissa de que os polos são conjugados.

# Exemplo

## Polinômio -> (s - 1)(s - 2) = s² - 3s + 2

```
raizes = [1.0, 2.0]
coefs  = coeficientesPolinomiais(raizes)
```

## Retorna -> [2.0, -3.0, 1.0] (2*s^0 - 3*s^1 + 1*s^2)
"""
function coeficientesPolinomiais(raizes::Vector{<:Number})
    n = length(raizes)

    # Inicializamos o vetor de coeficientes.
    poly      = zeros(ComplexF64, n + 1)
    poly_temp = zeros(ComplexF64, n + 1)

    poly[1] = 1.0 # Polinômio inicial: P(s) = 1 (ou seja, 1 * s^0) 

    for r in raizes
        poly_temp .= 0.0
        
        # Isso desloca todos os coeficientes um grau para cima (a_i * s se torna a_i * s^(i+1))
        for i in 1:(length(poly) - 1)
            poly_temp[i+1] += poly[i]
        end
        
        # Isso mantém o grau, mas escala o coeficiente por -r
        for i in eachindex(poly)
            poly_temp[i] -= r * poly[i]
        end
        
        poly .= poly_temp
    end
    return real.(poly)
end