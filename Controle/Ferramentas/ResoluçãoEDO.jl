using OrdinaryDiffEq
using DiffEqCallbacks

# =========================================================================
# CALLBACKS
    """
        callbackSolo(idx_Z=3, idx_Vz=6) -> ContinuousCallback

        Impede que a altitude (x[idx_Z]) fique negativa durante a simulação.
        Intercepta o cruzamento de Z = 0 e zera Z e Vz no instante exato,
        antes que qualquer passo seguinte seja calculado com Z < 0.

        Argumentos:
            - idx_Z  : índice do estado de altitude no vetor x (padrão: 3)
            - idx_Vz : índice da velocidade vertical no vetor x (padrão: 6)
    """
    function callbackSolo(; idx_Z=3, idx_Vz=6)
        condicao(x, t, integrador)  = x[idx_Z]

        function afeto!(integrador)
            integrador.u[idx_Z]  = max(0.0, integrador.u[idx_Z])
            integrador.u[idx_Vz] = max(0.0, integrador.u[idx_Vz])
        end

        return ContinuousCallback(condicao, afeto!,
                                rootfind    = true,
                                affect_neg! = afeto!)
    end

    """
        callbackHistoricoControle() -> SavingCallback

        Registra o sinal de controle u(t) a cada passo salvo pelo solver.
        Requer que o parâmetro do problema seja uma MalhaFechada com campos
        dim_planta, controlador e referencia.

        Retorna uma tupla (cb, historico_t, historico_u) para que os vetores
        de histórico possam ser acessados após a simulação.
    """
    function callbackHistoricoControle()
        historico_u = Vector{Vector{Float64}}()
        historico_t = Vector{Float64}()

        cb = SavingCallback(
            (x, t, integrador) -> begin
                n        = integrador.p.dim_planta
                x_planta = x[1:n]
                x_ctrl   = x[n+1:end]
                u        = calcularSaida(integrador.p.controlador,
                                        x_planta, x_ctrl,
                                        integrador.p.referencia, t)
                u_vec = u isa AbstractVector ? u : [u]
                push!(historico_t, t)
                push!(historico_u, copy(u_vec))
                nothing
            end,
            SavedValues(Float64, Nothing)
        )

        return cb, historico_t, historico_u
    end
# =========================================================================

# =========================================================================
# SOLVER
    """
        resolverSistema(sistema, x0, intervaloTempo, p; resolucao, salvar_controle, solo) -> solucao

        Resolve numericamente um sistema dinâmico no intervalo de tempo dado.

        Argumentos:
            - sistema          : função de dinâmica no formato (dx, x, p, t)
            - x0               : condições iniciais
            - intervaloTempo   : tupla (t_inicio, t_fim)
            - p                : parâmetros passados ao sistema (padrão: sistema)
            - resolucao        : passo de salvamento da solução (padrão: 0.01)
            - salvar_controle  : se true, registra u(t) ao longo da simulação
            - solo             : se true, ativa a restrição de altitude Z ≥ 0
    """
    function resolverSistema(sistema, x0, intervaloTempo, p=sistema; 
                            resolucao       = 0.01,
                            salvar_controle = false,
                            solo            = false)

        x0 = x0 isa Vector{Float64} ? x0 : collect(Float64, x0)

        # Monta lista de callbacks ativos
        callbacks = []

        historico_u = Vector{Vector{Float64}}()
        historico_t = Vector{Float64}()

        if salvar_controle
            cb_historico, historico_t, historico_u = callbackHistoricoControle()
            push!(callbacks, cb_historico)
        end

        if solo
            push!(callbacks, callbackSolo())
        end

        cb_final = if isempty(callbacks)
            nothing
        elseif length(callbacks) == 1
            callbacks[1]
        else
            CallbackSet(callbacks...)
        end

        problema = if p === nothing
            ODEProblem(sistema, x0, intervaloTempo)
        else
            ODEProblem(sistema, x0, intervaloTempo, p)
        end

        solucao = if cb_final === nothing
            solve(problema, saveat=resolucao)
        else
            solve(problema, saveat=resolucao, callback=cb_final)
        end

        if salvar_controle
            return (solucao=solucao, t_u=historico_t, u=historico_u)
        else
            return solucao
        end
    end
# =========================================================================