include("../Abstrações/Interfaces.jl")

# =========================================================================
# ESTRUTURA
    """
        Controlador PID clássico.

        Campos:
            - Kp : Ganho proporcional
            - Ki : Ganho integral
            - Kd : Ganho derivativo

        A referência e a planta são conectadas externamente via conectar().
    """
    struct PID <: Controlador
        Kp ::Float64
        Ki ::Float64
        Kd ::Float64

        function PID(; Kp, Ki, Kd)
            new(Kp, Ki, Kd)
        end
    end
# =========================================================================

# =========================================================================
# IMPLEMENTAÇÃO 
    numEstadosControle(::PID) = 1 # O PID tem 1 estado interno: o integrador do erro.

    #=
        x_planta[1] : saída medida 
        x_planta[2] : derivada da saída 
        x_ctrl[1]   : integral do erro acumulada
        ref         : vetor de referência — usa ref[1]
    =#
    function calcularSaida(pid::PID, x_planta, x_ctrl, ref, t)
        erro     = ref[1] - x_planta[1]
        integral = x_ctrl[1]
        derivada = -x_planta[2]        

        return pid.Kp * erro +
               pid.Ki * integral +
               pid.Kd * derivada
    end
    
    # O único estado interno é o integrador — sua derivada é simplesmente o erro.
    function evoluirEstado!(pid::PID, dx_ctrl, x_planta, x_ctrl, ref, t)
        dx_ctrl[1] = ref[1] - x_planta[1]  
    end
# =========================================================================