using ForwardDiff
using LinearAlgebra

"""
    Lineariza um sistema não-linear do FrameWork Controle ao redor 
    de um ponto de operação (x0, u0).

    # Argumentos:
    - `planta`: Objeto do tipo Planta (ex: PenduloInvertido) que contém os parâmetros e a função dinamica!.
    - `x0`: Vetor de estado no ponto de equilíbrio.
    - `u0`: Vetor de ação de controle no ponto de equilíbrio.
    - `h`: (Opcional) Função de saída `y = h(x, u)`. Se omitida, assume realimentação de todos os estados (y = x).

    # Retorno:
    - Matrizes `A`, `B`, `C` e `D` do Espaço de Estados (Linear).
"""
function linearizar(planta, x0::Vector{Float64}, u0::Vector{Float64}; h::Union{Function, Nothing}=nothing)
    
    # =======================================================================
    # O PULO DO GATO: Adaptador para o ForwardDiff
    # O ForwardDiff passa números especiais (Dual Numbers) para calcular a 
    # derivada exata. Nossa função precisa criar um vetor `dx` que suporte 
    # esses números, caso contrário ocorrerá um "Type Error".
    # =======================================================================
    function f_dinamica(x, u)
        # Descobre dinamicamente o tipo necessário (Float64 ou ForwardDiff.Dual)
        T = promote_type(eltype(x), eltype(u))
        
        # Cria o vetor dx com o tipo correto e o tamanho da planta
        dx = zeros(T, planta.numEstados)
        
        # O framework espera u escalar para SISO e vetor para MIMO.
        u_in = length(u) == 1 ? u[1] : u
        
        # Chama a dinâmica do framework (in-place)
        # t = 0.0 é arbitrário, pois os sistemas são invariantes no tempo (LTI)
        planta.dinamica!(dx, x, planta, 0.0; u=u_in)
        
        return dx
    end

    # Matriz A: Jacobiano da dinâmica em relação aos estados (mantendo u fixo em u0)
    A = ForwardDiff.jacobian(x -> f_dinamica(x, u0), x0)
    
    # Matriz B: Jacobiano da dinâmica em relação ao controle (mantendo x fixo em x0)
    B = ForwardDiff.jacobian(u -> f_dinamica(x0, u), u0)
    
    # Matrizes de Saída (C e D)
    if h === nothing
        n_estados = length(x0)
        n_entradas = length(u0)
        C = Matrix{Float64}(I, n_estados, n_estados)
        D = zeros(Float64, n_estados, n_entradas)
    else
        C = ForwardDiff.jacobian(x -> h(x, u0), x0)
        D = ForwardDiff.jacobian(u -> h(x0, u), u0)
    end
    
    return A, B, C, D
end