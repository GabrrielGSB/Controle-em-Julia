# Controle/Ferramentas/PropriedadesSistemas.jl

using LinearAlgebra

# =====================================================================
# FUNÇÕES AUXILIARES (Construção das matrizes P e Q)
    """
        _mat_control(A, B)

    Calcula e retorna a Matriz de Controlabilidade para o par (A, B).
    """
    function _mat_control(A::AbstractMatrix, B::AbstractMatrix)
        n = size(A, 1) 
        
        C_mat = copy(B)
        AB_k  = copy(B)
        
        # Concatena horizontalmente [B  AB  A^2B ... A^(n-1)B]
        for k in 1:(n-1)
            AB_k = A * AB_k
            C_mat = hcat(C_mat, AB_k)
        end
        
        return C_mat
    end

    """
        _mat_observ(A, C)

    Calcula e retorna a Matriz de Observabilidade para o par (A, C).
    """
    function _mat_observ(A::AbstractMatrix, C::AbstractMatrix)
        n = size(A, 1) 
        
        O_mat = copy(C)
        CA_k = copy(C)
        
        # Concatena verticalmente [C; CA; CA^2; ... ; CA^(n-1)]
        for k in 1:(n-1)
            CA_k = CA_k * A
            O_mat = vcat(O_mat, CA_k)
        end
        
        return O_mat
    end
# =====================================================================

# =====================================================================
# FUNÇÕES 
    """
        controlavel(A, B; tol=1e-7)

    Verifica se o sistema definido pelo par (A, B) é controlável.
    Retorna `true` se for controlável, `false` caso contrário.
    """
    function controlavel(A::AbstractMatrix, B::AbstractMatrix; retornar=false, tol=1e-7)
        n = size(A, 1)
        P = _mat_control(A, B)
        
        posto = rank(P, rtol=tol)

        if (retornar) return P
        else          return posto == n end
    end



    """
        observavel(A, C; retornar, tol=1e-7)

    Verifica se o sistema definido pelo par (A, C) é observável.
    Retorna `true` se for observável, `false` caso contrário.
    """
    function observavel(A::AbstractMatrix, C::AbstractMatrix; retornar=false, tol=1e-7)
        n = size(A, 1)
        Q = _mat_observ(A, C)
        
        posto = rank(Q, rtol=tol)
        
        if (retornar) return Q 
        else          return posto == n end
    end
# =====================================================================