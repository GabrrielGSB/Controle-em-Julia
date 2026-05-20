# =========================================================================
# GERAÇÃO DE POLOS POR ESPECIFICAÇÕES DE DESEMPENHO
# =========================================================================

"""
    projetarPolos(n::Int, Mp_percent::Float64, ts::Float64; fator_distancia::Float64 = 10.0)

Gera um vetor de `n` polos desejados a partir de especificações clássicas de desempenho 
da resposta transitória (aproximação de 2ª ordem).

Argumentos:
- `n`: Dimensão do sistema (número de variáveis de estado).
- `Mp_percent`: Sobressinal máximo desejado em percentagem (ex: 15.0 para 15%).
- `ts`: Tempo de acomodação desejado em segundos (usando o critério de 2%).
- `fator_distancia`: (Opcional) Quão mais rápidos os polos não-dominantes devem ser 
                     em relação à parte real dos polos dominantes (padrão: 10x).

Retorno:
- Um vetor de `ComplexF64` com `n` polos prontos para usar no Bass-Gura.
"""
function projetarPolos(n::Int, Mp_percent::Float64, ts::Float64; fator_distancia::Float64 = 10.0)
    # Validações de segurança
    if Mp_percent <= 0.0 || Mp_percent >= 100.0
        error("Erro: O sobressinal (Mp) deve estar estritamente entre 0 e 100.")
    end
    if ts <= 0.0
        error("Erro: O tempo de acomodação (ts) deve ser positivo.")
    end
    if n < 2
        error("Erro: O sistema deve ter pelo menos dimensão 2 para usar a aproximação de segunda ordem.")
    end

    # 1. Tradução Matemática: Desempenho Físico -> Parâmetros (zeta, wn)
    Mp = Mp_percent / 100.0
    zeta = -log(Mp) / sqrt(pi^2 + log(Mp)^2)
    wn = 4.0 / (zeta * ts)  # Fórmula clássica para tempo de acomodação de 2%

    # 2. Geração dos Polos Dominantes (Complexos Conjugados)
    sigma = -zeta * wn
    wd = wn * sqrt(1 - zeta^2)

    p1 = complex(sigma, wd)
    p2 = complex(sigma, -wd)

    # Inicializamos o vetor de polos com o par dominante
    polos = ComplexF64[p1, p2]

    # 3. Geração dos Polos Não-Dominantes (para n > 2, ex: Pêndulo Invertido)
    if n > 2
        # Definimos uma localização no semiplano esquerdo bem distante da origem
        sigma_longe = sigma * fator_distancia

        for i in 1:(n - 2)
            # Adicionamos uma pequeníssima separação entre os polos reais para evitar raízes 
            # exatamentes idênticas. Matrizes com autovalores repetidos podem causar 
            # singularidades numéricas em alguns solvers do Julia.
            polo_extra = complex(sigma_longe * (1.0 + 0.05 * (i - 1)), 0.0)
            push!(polos, polo_extra)
        end
    end

    return polos
end