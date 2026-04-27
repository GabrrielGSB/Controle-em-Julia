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
# IMPLEMENTAÇÃO SISO
    numEstadosControle(::PID) = 1

    function calcularSaida(pid::PID, x_planta, x_ctrl, ref, t)
        erro     = ref[1] - x_planta[1]
        integral = x_ctrl[1]
        derivada = -x_planta[2]        

        return pid.Kp * erro +
               pid.Ki * integral +
               pid.Kd * derivada
    end
    
    function evoluirEstado!(pid::PID, dx_ctrl, x_planta, x_ctrl, ref, t)
        dx_ctrl[1] = ref[1] - x_planta[1]  
    end
# =========================================================================

# =========================================================================
# IMPLEMENTAÇÃO MIMO 
    """
        PID_MIMO

        Orquestrador de N controladores PID independentes para sistemas MIMO.
        Cada canal tem seu próprio PID e integrador.

        Campos:
            - controladores : Vetor de PID, um por entrada
            - ix_saida      : Índices dos estados de saída medida em x_planta  (ex: [3, 7, 8, 9])
            - idx_saida     : Índices dos estados de derivada da saída em x_planta (ex: [6, 10, 11, 12])
            - mapear_saida  : Função opcional (u_canais, x_planta) -> u_planta para transformar
                              as saídas dos PIDs antes de entregar à planta (ex: feedforward de gravidade)

        Exemplo — 4 PIDs para drone:
            PID_MIMO(
                controladores = [pid_Z, pid_roll, pid_pitch, pid_yaw],
                ix_saida      = [3, 7, 8, 9],
                idx_saida     = [6, 10, 11, 12],
                mapear_saida  = (u, x) -> [m*g + u[1], u[2], u[3], u[4]]
            )
    """
    struct PID_MIMO <: Controlador
        controladores ::Vector{PID}
        ix_saida      ::Vector{Int}
        idx_saida     ::Vector{Int}
        mapear_saida  ::Function

        function PID_MIMO(; controladores, ix_saida, idx_saida, mapear_saida = (u, x) -> u)
            @assert length(controladores) == length(ix_saida) == length(idx_saida) """
            Número de controladores ($(length(controladores))), índices de saída ($(length(ix_saida))) 
            e índices de derivada ($(length(idx_saida))) devem ser iguais.
            """

            new(controladores, ix_saida, idx_saida, mapear_saida)
        end
    end

    # Um integrador por canal
    numEstadosControle(ctrl::PID_MIMO) = length(ctrl.controladores)

    function calcularSaida(ctrl::PID_MIMO, x_planta, x_ctrl, ref, t)
        n = length(ctrl.controladores)
        u = Vector{Float64}(undef, n)

        for i in 1:n
            x_sys = [x_planta[ctrl.ix_saida[i]], x_planta[ctrl.idx_saida[i]]]
            u[i] = calcularSaida(ctrl.controladores[i], x_sys, [x_ctrl[i]], [ref[i]], t)
        end

        return ctrl.mapear_saida(u, x_planta)
    end

    function evoluirEstado!(ctrl::PID_MIMO, dx_ctrl, x_planta, x_ctrl, ref, t)
        for i in 1:length(ctrl.controladores)
            dx_ctrl[i] = ref[i] - x_planta[ctrl.ix_saida[i]]
        end
    end
# =========================================================================