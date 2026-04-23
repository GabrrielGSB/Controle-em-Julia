# Controle/Ferramentas/AnálisePerformance.jl
#
# Módulo de análise de performance de controladores.
# Compatível com qualquer solução retornada por resolverSistema().
#
# USO BÁSICO:
#   metricas = analisarPerformance(solucao, referencia=1.0, idx_estado=3, idx_controle=nothing)
#   imprimirRelatorio(metricas)
#
# USO AVANÇADO (com sinal de controle explícito):
#   metricas = analisarPerformance(solucao, u_historico, referencia=1.0, idx_estado=3)

using Statistics
using Printf

# =========================================================================
# STRUCT DE RESULTADOS
# =========================================================================

"""
    Metricas

    Contém todas as métricas de performance calculadas para um controlador.
    Campos organizados em grupos temáticos:

        RASTREAMENTO:
            - tempo_subida       : Tempo de 10% a 90% da referência [s]
            - tempo_pico         : Tempo até o valor de pico [s]
            - tempo_assentamento : Tempo até ficar em ±banda% da ref [s]
            - sobressinal_pct    : Pico máximo acima da referência [%]
            - erro_regime        : Erro residual em regime permanente
            - banda_assentamento : Banda usada no cálculo (padrão 0.02 = 2%)

        ÍNDICES INTEGRAIS (calculados sobre o erro e(t) = ref - y(t)):
            - ISE   : Integral do Quadrado do Erro  ∫ e² dt
            - IAE   : Integral do Valor Absoluto    ∫ |e| dt
            - ITAE  : Integral Tempo × |e|          ∫ t·|e| dt
            - ITSE  : Integral Tempo × e²           ∫ t·e² dt

        ESFORÇO DE CONTROLE (requer u_historico):
            - energia_controle   : Norma L2 do controle — ∫ u² dt
            - variacao_total     : TV do controle — ∫ |du/dt| dt
            - u_max              : Valor máximo absoluto do controle
            - u_rms              : Valor RMS do sinal de controle

        ROBUSTEZ / TRANSITÓRIO:
            - num_oscilacoes     : Número de cruzamentos do sinal pelo valor de referência
            - decaimento         : Taxa de decaimento exponencial estimada
"""
struct Metricas
    # --- Rastreamento ---
    tempo_subida        ::Union{Float64, Nothing}
    tempo_pico          ::Union{Float64, Nothing}
    tempo_assentamento  ::Union{Float64, Nothing}
    sobressinal_pct     ::Float64
    erro_regime         ::Float64
    banda_assentamento  ::Float64

    # --- Índices Integrais ---
    ISE    ::Float64
    IAE    ::Float64
    ITAE   ::Float64
    ITSE   ::Float64

    # --- Esforço de Controle ---
    energia_controle    ::Union{Float64, Nothing}
    variacao_total      ::Union{Float64, Nothing}
    u_max               ::Union{Float64, Nothing}
    u_rms               ::Union{Float64, Nothing}

    # --- Robustez / Transitório ---
    num_oscilacoes      ::Int
    decaimento          ::Union{Float64, Nothing}

    # --- Metadados ---
    referencia          ::Float64
    idx_estado          ::Int
    duracao_total       ::Float64
end

# =========================================================================
# FUNÇÃO PRINCIPAL

"""
    analisarPerformance(solucao; referencia, idx_estado, u_historico, banda) -> Metricas

    Calcula todas as métricas de performance para um sistema simulado.

    Argumentos:
        - solucao     : Objeto de solução retornado por resolverSistema()
        - referencia  : Valor de referência (setpoint) desejado
        - idx_estado  : Índice do estado a ser avaliado (padrão: 1)
        - u_historico : (opcional) Vetor Float64 do sinal de controle u(t) ao longo do tempo
                        Deve ter o mesmo comprimento que solucao.t
        - banda       : Tolerância para tempo_assentamento (padrão: 0.02 = ±2%)

    Exemplo:
        sol = resolverSistema(pendulo_malha_fechada!, x0, (0.0, 10.0), params)
        m   = analisarPerformance(sol, referencia=π, idx_estado=1)
        imprimirRelatorio(m)
"""
function analisarPerformance(solucao;
                              referencia  ::Real    = 0.0,
                              idx_estado  ::Int     = 1,
                              u_historico ::Union{Vector{Float64}, Nothing} = nothing,
                              banda       ::Float64 = 0.02)

    t    = solucao.t
    y    = [u[idx_estado] for u in solucao.u]
    ref  = Float64(referencia)
    dt   = diff(t)  
    n    = length(t)
    x0   = y[1]
    span = ref - x0

    # === ERRO ===
    e_t     = ref .- y     
    e_t_abs = abs.(e_t)
    e_t_sq  = e_t .^ 2

    # === ÍNDICES INTEGRAIS ===
    ISE  = _integrar(e_t_sq,       t)
    IAE  = _integrar(e_t_abs,      t)
    ITAE = _integrar(t .* e_t_abs, t)
    ITSE = _integrar(t .* e_t_sq,  t)

    # === MÉTRICAS DE RASTREAMENTO ===
    tempo_assentamento = _calcularTempoAssentamento(t, y, ref, banda)
    sobressinal_pct    = _calcularSobressinal(y, ref)
    num_oscilacoes     = _contarOscilacoes(y, ref)
    tempo_subida       = _calcularTempoSubida(t, y, x0, ref, span)
    erro_regime        = _calcularErroRegime(y, ref)
    tempo_pico         = _calcularTempoPico(t, y, span)
    decaimento         = _estimarDecaimento(t, e_t_abs)

    # === ESFORÇO DE CONTROLE ===
    energia_controle = nothing
    variacao_total   = nothing
    u_max            = nothing
    u_rms            = nothing

    if u_historico !== nothing
        @assert length(u_historico) == n "u_historico deve ter o mesmo comprimento que solucao.t"

        u = u_historico
        du = diff(u) ./ dt

        energia_controle = _integrar(u .^ 2, t)
        variacao_total   = _integrar(abs.(du), t[1:end-1])
        u_max            = maximum(abs.(u))
        u_rms            = sqrt(_integrar(u .^ 2, t) / t[end])
    end

    return Metricas(
        tempo_subida,
        tempo_pico,
        tempo_assentamento,
        sobressinal_pct,
        erro_regime,
        banda,
        ISE, IAE, ITAE, ITSE,
        energia_controle,
        variacao_total,
        u_max,
        u_rms,
        num_oscilacoes,
        decaimento,
        ref,
        idx_estado,
        t[end]
    )
end

# =========================================================================
# COMPARAÇÃO DE MÚLTIPLOS CONTROLADORES
# =========================================================================

"""
    compararControladores(solucoes, nomes; referencia, idx_estado, ...) -> Vector{Metricas}

    Analisa e compara múltiplos controladores lado a lado.

    Exemplo:
        m_pid = analisarPerformance(sol_pid, referencia=1.0)
        m_lqr = analisarPerformance(sol_lqr, referencia=1.0)
        compararControladores([sol_pid, sol_lqr], ["PID", "LQR"], referencia=1.0)
"""
function compararControladores(solucoes::AbstractVector,
                                nomes::Vector{String};
                                kwargs...)
    metricas = [analisarPerformance(sol; kwargs...) for sol in solucoes]
    imprimirComparacao(metricas, nomes)
    return metricas
end

# =========================================================================
# IMPRESSÃO NO TERMINAL
# =========================================================================

"""
    imprimirRelatorio(m::Metricas; nome)

    Imprime um relatório formatado no terminal com todas as métricas disponíveis.
"""
function imprimirRelatorio(m::Metricas; nome::String = "Controlador")
    linha = "="^56
    sublinha = "-"^56
    _fmt(v::Float64) = @sprintf("%.6g", v)
    _fmt(::Nothing) = "N/D (u_historico não fornecido)"

    println("\n", linha)
    println("  RELATÓRIO DE PERFORMANCE — $nome")
    println("  Estado analisado : y[$(m.idx_estado)]")
    println("  Referência       : $(m.referencia)")
    println("  Duração          : $(m.duracao_total) s")
    println(sublinha)

    println("  RASTREAMENTO DE REFERÊNCIA")
    println("    Tempo de subida       : $(_fmt_t(m.tempo_subida)) s")
    println("    Tempo de pico         : $(_fmt_t(m.tempo_pico)) s")
    println("    Tempo assentamento    : $(_fmt_t(m.tempo_assentamento)) s   (banda ±$(round(Int, m.banda_assentamento*100))%)")
    println("    Sobressinal           : $(round(m.sobressinal_pct, digits=3)) %")
    println("    Erro regime perm.     : $(_fmt(m.erro_regime))")
    println(sublinha)

    println("  ÍNDICES INTEGRAIS DE ERRO")
    println("    ISE  (∫ e²·dt)        : $(_fmt(m.ISE))")
    println("    IAE  (∫|e|·dt)        : $(_fmt(m.IAE))")
    println("    ITAE (∫t·|e|·dt)      : $(_fmt(m.ITAE))")
    println("    ITSE (∫t·e²·dt)       : $(_fmt(m.ITSE))")
    println(sublinha)

    println("  ESFORÇO DE CONTROLE")
    println("    Energia (∫ u²·dt)     : $(_fmt(m.energia_controle))")
    println("    Variação total (TV)   : $(_fmt(m.variacao_total))")
    println("    |u| máximo            : $(_fmt(m.u_max))")
    println("    u RMS                 : $(_fmt(m.u_rms))")
    println(sublinha)

    println("  QUALIDADE DO TRANSITÓRIO")
    println("    Nº de oscilações      : $(m.num_oscilacoes)")
    println("    Taxa de decaimento    : $(_fmt_t(m.decaimento))")
    println(linha, "\n")
end

"""
    Imprime uma tabela comparativa de múltiplos controladores.
"""
function imprimirComparacao(metricas::Vector{Metricas}, nomes::Vector{String})
    println("\n", "="^80)
    println("  COMPARAÇÃO DE CONTROLADORES")
    println("="^80)

    col = 15
    header = rpad("Métrica", 28)
    for n in nomes
        header *= rpad(n, col)
    end
    println(header)
    println("-"^80)

    linhas = [
        ("Tempo subida (s)",     m -> _fmt_t(m.tempo_subida)),
        ("Tempo assent. (s)",    m -> _fmt_t(m.tempo_assentamento)),
        ("Sobressinal (%)",      m -> @sprintf("%.3g", m.sobressinal_pct)),
        ("Erro regime",          m -> @sprintf("%.4g", m.erro_regime)),
        ("ISE",                  m -> @sprintf("%.4g", m.ISE)),
        ("IAE",                  m -> @sprintf("%.4g", m.IAE)),
        ("ITAE",                 m -> @sprintf("%.4g", m.ITAE)),
        ("ITSE",                 m -> @sprintf("%.4g", m.ITSE)),
        ("Energia controle",     m -> _fmt(m.energia_controle)),
        ("Variação total (TV)",  m -> _fmt(m.variacao_total)),
        ("|u| máximo",           m -> _fmt(m.u_max)),
        ("Nº oscilações",        m -> string(m.num_oscilacoes)),
        ("Decaimento",           m -> _fmt_t(m.decaimento)),
    ]

    for (label, fn) in linhas
        linha = rpad(label, 28)
        for m in metricas
            linha *= rpad(fn(m), col)
        end
        println(linha)
    end
    println("="^80, "\n")
end

# =========================================================================
# EXPORTAR PARA DICT (útil para plotagem ou logging)
# =========================================================================

"""
    paraDict(m::Metricas) -> Dict{String, Any}

Converte as métricas para um dicionário, útil para plotagem ou exportação.
"""
function paraDict(m::Metricas)
    return Dict(
        "tempo_subida"       => m.tempo_subida,
        "tempo_pico"         => m.tempo_pico,
        "tempo_assentamento" => m.tempo_assentamento,
        "sobressinal_pct"    => m.sobressinal_pct,
        "erro_regime"        => m.erro_regime,
        "ISE"                => m.ISE,
        "IAE"                => m.IAE,
        "ITAE"               => m.ITAE,
        "ITSE"               => m.ITSE,
        "energia_controle"   => m.energia_controle,
        "variacao_total"     => m.variacao_total,
        "u_max"              => m.u_max,
        "u_rms"              => m.u_rms,
        "num_oscilacoes"     => m.num_oscilacoes,
        "decaimento"         => m.decaimento,
    )
end

# =========================================================================
# FUNÇÕES AUXILIARES INTERNAS

    function _calcularErroRegime(y::Vector{Float64}, ref::Float64)
        n = length(y)
        
        janela_final = y[max(1, end-round(Int, 0.1 * n)):end]
        y_final = mean(janela_final)
        
        return ref - y_final
    end

    # Porcentagem do sobressinal em relação a referência
    function _calcularSobressinal(y::Vector{Float64}, ref::Float64)
        if ref == 0.0
            return 100.0 * maximum(abs.(y))
        
        elseif ref > 0.0
            pico = maximum(y)
            if pico > ref
                return 100.0 * (pico - ref) / ref
            else
                return 0.0
            end
            
        else
            pico_negativo = minimum(y)
            if pico_negativo < ref
                return 100.0 * (ref - pico_negativo) / abs(ref)
            else
                return 0.0
            end
        end
    end

    # Integração numérica pela regra dos trapézios
    function _integrar(y::Vector{Float64}, t::Vector{Float64})
        n = min(length(y), length(t))
        total = 0.0
        for i in 1:(n-1)
            total += 0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i])
        end
        return total
    end

    # Tempo de subida: de 10% a 90% do span (x0 → ref)
    function _calcularTempoSubida(t, y, x0, ref, span)
        abs(span) < 1e-10 && return nothing
        thr_low  = x0 + 0.10 * span
        thr_high = x0 + 0.90 * span
        t_low = t_high = nothing
        for i in 1:length(t)
            if t_low === nothing && sign(span) * y[i] >= sign(span) * thr_low
                t_low = t[i]
            end
            if t_high === nothing && sign(span) * y[i] >= sign(span) * thr_high
                t_high = t[i]
                break
            end
        end
        (t_low === nothing || t_high === nothing) && return nothing
        return t_high - t_low
    end

    function _calcularTempoPico(t::Vector{Float64}, y::Vector{Float64}, span::Float64)
        idx_pico = argmax(sign(span) .* y)
        return t[idx_pico]
    end

    # Tempo de assentamento: último instante em que |e| > banda * |ref|
    function _calcularTempoAssentamento(t, y, ref, banda)
        tolerancia = banda * abs(ref)
        tolerancia < 1e-10 && (tolerancia = banda)  # ref = 0 fallback
        ultimo_fora = nothing
        for i in 1:length(t)
            if abs(y[i] - ref) > tolerancia
                ultimo_fora = i
            end
        end
        ultimo_fora === nothing && return t[1]
        ultimo_fora == length(t) && return nothing   # Nunca assentou
        return t[ultimo_fora]
    end

    # Conta cruzamentos pela referência (cada par cruza = 1 oscilação)
    function _contarOscilacoes(y, ref)
        sinais = sign.(y .- ref)
        count = 0
        for i in 2:length(sinais)
            if sinais[i] != sinais[i-1] && sinais[i] != 0
                count += 1
            end
        end
        return count ÷ 2
    end

    # Estima taxa de decaimento via regressão linear no log do envelope
    function _estimarDecaimento(t, e_t_abs)
        eps = 1e-12
        log_e = log.(e_t_abs .+ eps)
        # Regressão linear: log(e) ≈ a + b*t  →  b é a taxa de decaimento
        n = length(t)
        t_mean = mean(t)
        e_mean = mean(log_e)
        num = sum((t .- t_mean) .* (log_e .- e_mean))
        den = sum((t .- t_mean) .^ 2)
        abs(den) < 1e-10 && return nothing
        b = num / den
        return b   # negativo = sistema decaindo (estável)
    end
# =========================================================================

_fmt(::Nothing)  = "N/D"
_fmt(v::Float64) = @sprintf("%.4g", v)
_fmt_t(::Nothing)  = "N/D"
_fmt_t(v::Float64) = @sprintf("%.4g", v)