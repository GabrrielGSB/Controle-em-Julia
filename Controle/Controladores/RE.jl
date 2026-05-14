include("../Abstrações/Interfaces.jl")

# =========================================================================
# STRUCT DO CONTROLADOR
# =========================================================================

"""
    Controlador de Realimentação de Estados Completo (State Feedback).
    
    Campos:
        - K : Matriz de ganhos de realimentação (AbstractMatrix para suportar adjuntas/transpostas).
              A lei de controle assume a forma u = -K * (x - ref).
"""
struct RealimentacaoEstado <: Controlador
    K::AbstractMatrix{Float64}
end

# =========================================================================
# INTERFACE
# =========================================================================

"""
    dimEstado(::RealimentacaoEstado)
    A realimentação de estados pura é um controlador estático, 
    não possui estados internos (memória).
"""
numEstadosControle(::RealimentacaoEstado) = 0

"""
    calcularSaida(ctrl, x_planta, x_ctrl, ref, t) -> u
    
    Aplica a lei de controle u = -K * (x_planta - ref).
    Se o objetivo for apenas regulação para a origem, 'ref' pode ser um vetor de zeros.
"""
function calcularSaida(ctrl::RealimentacaoEstado, x_planta, x_ctrl, ref, t)
    u_vec = -ctrl.K * (x_planta - ref)
    
    if (length(u_vec) == 1) return u_vec[1]
    else                    return u_vec
    end
end

"""
    evoluirEstado(ctrl, dx_ctrl, x_planta, x_ctrl, ref, t)
    
    Como dimEstado é 0, este controlador não tem derivadas internas para atualizar.
    A função não faz nada.
"""
function evoluirEstado(ctrl::RealimentacaoEstado, dx_ctrl, x_planta, x_ctrl, ref, t)
    return nothing
end