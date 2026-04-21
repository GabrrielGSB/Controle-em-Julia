# Controle/Abstrações/Interfaces.jl

# =========================================================================
# TIPOS ABSTRATOS
    """
        Supertipo de toda planta dinâmica do framework.
        Todo struct que represente um sistema físico deve herdar deste tipo.

        Exemplo:
            struct PenduloParams <: SistemaPlanta ... end
    """
    abstract type SistemaPlanta end

    """
        Supertipo de todo controlador do framework.
        Todo struct que represente uma lei de controle deve herdar deste tipo.

        Exemplo:
            struct PID <: SistemaControlador ... end
    """
    abstract type SistemaControlador end
# =========================================================================

# =========================================================================
# INTERFACE DO CONTROLADOR — funções que todo controlador deve implementar
    """
        calcularSaida(controlador, x_planta, x_ctrl, ref, t) -> u

        Calcula a ação de controle u no instante t.

        Argumentos:
            - controlador : struct do controlador (PID, LQR, MPC...)
            - x_planta    : vetor de estados atual da planta
            - x_ctrl      : vetor de estados internos do controlador (ex: integrador)
            - ref         : referência / setpoint (vindo do usuário ou de outro controlador)
            - t           : tempo atual
    """
    function calcularSaida end

    """
        evoluirEstado(controlador, dx_ctrl, x_planta, x_ctrl, ref, t)

        Escreve em `dx_ctrl` a derivada dos estados internos do controlador.
        Para controladores sem memória (LQR, SOS, Fuzzy), deve ser definido
        mas pode ter corpo vazio.

        Argumentos: mesmos de calcularSaida, mais dx_ctrl (vetor de saída, in-place).
    """
    function evoluirEstado end

    """
        dimEstado(controlador) -> Int

        Retorna a dimensão do vetor de estados internos do controlador.
            - PID       → 1  (o integrador do erro)
            - LQR, SOS  → 0  (lei estática, sem memória)
            - LQG       → n  (dimensão do estado estimado)
    """
    function dimEstado end
# =========================================================================