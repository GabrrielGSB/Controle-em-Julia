# Controle/Abstrações/MalhaFechada.jl

include("Interfaces.jl")

# =========================================================================
# STRUCT DA MALHA FECHADA
    """
        Representa um sistema em malha fechada: planta + controlador conectados.

        O vetor de estado aumentado é organizado como:
        x = [x_planta (dim_planta); x_ctrl (dim_ctrl)]

        Criado via conectar() — não instanciar diretamente.
    """
    struct MalhaFechada{P <: SistemaPlanta, C <: SistemaControlador}
        planta      ::P
        controlador ::C
        referencia  ::Vector{Float64}
        dim_planta  ::Int
        dim_ctrl    ::Int
    end
# =========================================================================

# =========================================================================
# CONSTRUTOR — conectar()
    """
        conectar(planta, controlador, referencia) -> MalhaFechada

        Conecta uma planta a um controlador em malha fechada.
        Equivalente ao feedback() do MATLAB.

        Argumentos:
            - planta        : struct <: SistemaPlanta  (ex: PenduloParams)
            - controlador   : struct <: SistemaControlador (ex: PID)
            - referencia    : valor desejado — escalar ou vetor

        Exemplo:
            sys = conectar(pendulo, pid, π)
            sys = conectar(pendulo, lqr, [1.0, 0.0])
    """
    function conectar(planta::P, controlador::C, referencia) where {P <: SistemaPlanta, C <: SistemaControlador}
        ref = referencia isa Vector ? referencia : [Float64(referencia)]

        return MalhaFechada(planta,
                            controlador,
                            ref,
                            planta.dim,
                            dimEstado(controlador))
    end
# =========================================================================

# =========================================================================
# CONDIÇÕES INICIAIS — monta o vetor aumentado
    """
        condicoesIniciais(malha, x0_planta) -> Vector

        Monta o vetor de condições iniciais aumentado [x_planta; zeros(dim_ctrl)].
        Os estados internos do controlador começam sempre em zero.

        Exemplo:
        x0 = condicoesIniciais(sys, [1.0, 0.0])  # → [1.0, 0.0, 0.0] para PID
    """
    function condicoesIniciais(malha::MalhaFechada, x0_planta::Vector)
        @assert length(x0_planta) == malha.dim_planta """
        Dimensão incorreta: x0_planta tem $(length(x0_planta)) estados,
        mas a planta espera $(malha.dim_planta).
        """
        return [x0_planta; zeros(malha.dim_ctrl)]
    end
# =========================================================================

# =========================================================================
# DINÂMICA AUMENTADA — callable struct, compatível com ODEProblem
    """
        Dinâmica aumentada da malha fechada.
        Chamada pelo solver como malha(dx, x, p, t).

        Organização interna:
        x  = [x_planta | x_ctrl]
        dx = [dx_planta | dx_ctrl]
    """
    function (malha::MalhaFechada)(dx, x, _, t)
        n = malha.dim_planta
        
        # Fatia o vetor de estado sem alocação
        x_planta  = @view x[1:n]
        x_ctrl    = @view x[n+1:end]
        dx_planta = @view dx[1:n]
        dx_ctrl   = @view dx[n+1:end]
        
        # 1. Controlador lê a planta e calcula u
        u = calcularSaida(malha.controlador, x_planta, x_ctrl, malha.referencia, t)
        
        # 2. Planta evolui com u
        malha.planta.dinamica!(dx_planta, x_planta, malha.planta, t; u=u)
        
        # 3. Estados internos do controlador evoluem
        if malha.dim_ctrl > 0
            evoluirEstado(malha.controlador, dx_ctrl, x_planta, x_ctrl, malha.referencia, t)
        end
    end
# =========================================================================