using PlotlyJS 

"""
    Plota estados e/ou sinal de controle ao longo do tempo.

    Argumentos:
        - resultado        : Solução do ODE *ou* NamedTuple com campos (solucao, t_u, u)
                             retornado por resolverSistema(salvar_controle=true)
        - titulo           : Título do gráfico
        - estados          : Índices dos estados a plotar. `nothing` oculta os estados.
        - mostrar_controle : Se true, adiciona subplots do sinal de controle u(t).
                             Ignorado se resultado não contiver histórico de controle.
"""
function plotarNoTempo(resultado;
                       titulo           = "Análise Temporal do Sistema",
                       estados          = nothing,
                       mostrar_controle = false)

    # ── 1. Desempacota resultado (NamedTuple ou solução direta) ──────────────
    tem_controle = resultado isa NamedTuple         &&
                   hasproperty(resultado, :solucao) &&
                   hasproperty(resultado, :t_u)     &&
                   hasproperty(resultado, :u)

    solucao    = tem_controle ? resultado.solucao : resultado
    t_estados  = solucao.t

    # Estados padrão: todos, caso não especificado e não foram ocultados
    indices_estados = estados === nothing ? (1:length(solucao.u[1])) :
                      estados == 0        ? Int[]                    : 
                     (estados isa Int ? [estados] : collect(estados))

    # ── 2. Valida pedido de controle ─────────────────────────────────────────
    if mostrar_controle && !tem_controle
        error("Histórico de controle não encontrado. " *
              "Rode resolverSistema com salvar_controle=true.")
    end

    # ── 3. Monta lista de subplots ────────────────────────────────────────────
    #   Cada entrada: (tipo, índice, título)
    #   tipo = :estado ou :controle
    subplots = Tuple{Symbol, Int, String}[]

    for i in indices_estados
        push!(subplots, (:estado, i, "x$i"))
    end

    if mostrar_controle && tem_controle
        dim_u = length(resultado.u[1])
        for i in 1:dim_u
            push!(subplots, (:controle, i, "u$i"))
        end
    end

    num_graficos = length(subplots)
    num_graficos == 0 && error("Nenhum estado ou controle selecionado para plotar.")

    # ── 4. Cria a figura ──────────────────────────────────────────────────────
    titulos_subplots = reshape([sp[3] for sp in subplots], 1, num_graficos)

    fig = make_subplots(
        rows              = num_graficos,
        cols              = 1,
        shared_xaxes      = true,
        vertical_spacing  = 0.15 / max(1, num_graficos - 1),
        subplot_titles    = titulos_subplots
    )

    # ── 5. Preenche cada subplot ──────────────────────────────────────────────
    for (linha, (tipo, idx, nome)) in enumerate(subplots)
        if tipo == :estado
            dados = [u[idx] for u in solucao.u]
            t_plot = t_estados
            cor = attr(width=2)
        else  # :controle
            dados  = [vetor[idx] for vetor in resultado.u]
            t_plot = resultado.t_u
            cor = attr(width=2, dash="dot")
        end

        add_trace!(fig,
            scatter(x=t_plot, y=dados, mode="lines", name=nome, line=cor),
            row=linha, col=1)
    end

    # ── 6. Layout ─────────────────────────────────────────────────────────────
    config = Dict(
        :title_text  => titulo,
        :title_x     => 0.5,
        :height      => max(400, 250 * num_graficos),
        :width       => 800,
        :showlegend  => true,
        :hovermode   => "x unified"
    )
    config[Symbol("xaxis$(num_graficos)_title")] = "Tempo (s)"

    relayout!(fig; config...)
    display(fig)
end

"""
    Compara N sistemas no tempo. 
"""
function plotarNoTempo(solucoes::AbstractVector; 
                       nomes=nothing,
                       titulo = "Comparação Temporal de Sistemas", 
                       estados = 1:length(solucoes[1].u[1])) 
    
    num_sistemas = length(solucoes)
    num_graficos = length(estados)
    
    # 1. Ajuste automático de nomes se o usuário não passar a lista
    if nomes === nothing
        nomes = ["Sistema $i" for i in 1:num_sistemas]
    end
    
    # Paleta de cores para manter o mesmo sistema com a mesma cor em todos os gráficos
    paleta = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
    
    # 2. Cria os títulos dos subplots dinamicamente 
    titulos_subplots = reshape(["Evolução de x$i" for i in estados], 1, num_graficos)
    
    # 3. Cria o layout de subplots dinâmico
    fig = make_subplots(
        rows=num_graficos, cols=1, 
        shared_xaxes=true,      
        vertical_spacing=0.15 / max(1, (num_graficos - 1)),
        subplot_titles=titulos_subplots
    )
    
    # 4. Loops Aninhados: Para cada estado, plota a linha de cada sistema
    for (linha_atual, estado_idx) in enumerate(estados)
        
        for (idx_sol, sol) in enumerate(solucoes)
            t = sol.t
            x_dados = [u[estado_idx] for u in sol.u]
            
            trace = scatter(
                x=t, 
                y=x_dados, 
                mode="lines", 
                # O nome na legenda será "Sistema 1", "Sistema 2", etc.
                name="$(nomes[idx_sol])", 
                
                # A MÁGICA 1: legendgroup conecta as linhas. Se clicar na legenda, desliga todas as variáveis daquele sistema
                legendgroup="grupo_$(idx_sol)", 
                
                # A MÁGICA 2: Mostra o nome na legenda APENAS no primeiro subplot para não duplicar nomes
                showlegend=(linha_atual == 1),
                
                # Mantém a mesma cor para o mesmo sistema em todos os subplots
                line=attr(width=2, color=paleta[mod1(idx_sol, length(paleta))])
            )
            add_trace!(fig, trace, row=linha_atual, col=1)
        end
    end
    
    # 5. Configuração dinâmica do layout
    config_layout = Dict(
        :title_text => titulo, 
        :title_x => 0.5, 
        :height => max(400, 250 * num_graficos), 
        :width => 800,
        :hovermode => "x unified" # MÁGICA 3: O mouse mostra os valores de todos os sistemas no mesmo instante t
    )
    
    chave_eixo_x = Symbol("xaxis$(num_graficos)_title")
    config_layout[chave_eixo_x] = "Tempo (s)"
    
    relayout!(fig; config_layout...)
    display(fig)
end

"""
    Plota o retrato de fase geométrico entre duas variáveis de estado à sua escolha.
    - estados: Uma tupla com os índices das duas variáveis. Ex: (1, 3) plota x1 no eixo X e x3 no eixo Y.
"""
function plotarRetratoFase(solucao; 
                           estados=(1, 2), 
                           titulo="Retrato de Fase")
    
    # 1. Trava de segurança para o usuário
    if length(estados) != 2
        error("Um retrato de fase 2D precisa de exatamente 2 estados. Passe uma tupla como estados=(1, 3)")
    end
    
    idx_x, idx_y = estados
    
    # 2. Extração cirúrgica dos dados das variáveis escolhidas
    x_dados = [u[idx_x] for u in solucao.u]
    y_dados = [u[idx_y] for u in solucao.u]
    
    # 3. Criação das camadas do gráfico (Traces)
    
    # A. A trajetória do sistema
    trace_traj = scatter(
        x=x_dados, y=y_dados, 
        mode="lines", 
        name="Trajetória", 
        line=attr(width=2, color="#1f77b4") 
    )
    
    # B. Ponto de Partida (Condição Inicial)
    trace_inicio = scatter(
        x=[x_dados[1]], y=[y_dados[1]], 
        mode="markers", 
        name="Início", 
        marker=attr(size=10, color="green")
    )
    
    # C. O Alvo de Controle (A Origem 0,0)
    trace_origem = scatter(
        x=[0.0], y=[0.0], 
        mode="markers", 
        name="Origem", 
        marker=attr(size=12, color="black", symbol="x")
    )
    
    # 4. Configuração Geométrica do Layout
    layout = Layout(
        title_text=titulo,
        title_x=0.5,
        xaxis_title="Estado x$(idx_x)",
        yaxis_title="Estado x$(idx_y)",
        width=700,
        height=700,
        showlegend=true,
        
        yaxis_scaleanchor="x",
        yaxis_scaleratio=1
    )
    
    # 5. Junta tudo e exibe
    fig = plot([trace_traj, trace_inicio, trace_origem], layout)
    display(fig)
end

"""
    Compara o retrato de fase de múltiplos sistemas no mesmo plano geométrico.
"""
function plotarRetratoFase(solucoes::AbstractVector; 
                           estados=(1, 2), 
                           nomes=nothing,
                           titulo="Comparação de Retratos de Fase")
    
    # 1. Trava de segurança
    if length(estados) != 2
        error("Um retrato de fase 2D precisa de exatamente 2 estados. Ex: estados=(1, 2)")
    end
    
    idx_x, idx_y = estados
    num_sistemas = length(solucoes)
    
    # Nomes automáticos se não fornecidos
    if nomes === nothing
        nomes = ["Sistema $i" for i in 1:num_sistemas]
    end
    
    # Paleta de cores oficial
    paleta = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
    
    # Vetor vazio para guardar todas as "camadas" do gráfico
    traces = GenericTrace[] 
    
    # 2. Loop principal: Cria as linhas e pontos para cada sistema
    for (idx_sol, sol) in enumerate(solucoes)
        cor_sistema = paleta[mod1(idx_sol, length(paleta))]
        nome_sistema = nomes[idx_sol]
        
        x_dados = [u[idx_x] for u in sol.u]
        y_dados = [u[idx_y] for u in sol.u]
        
        # A. A trajetória do sistema
        push!(traces, scatter(
            x=x_dados, y=y_dados, 
            mode="lines", 
            name=nome_sistema, 
            legendgroup="grupo_$idx_sol", # MÁGICA: Conecta a linha com os pontos na legenda
            line=attr(width=2, color=cor_sistema)
        ))
        
        # B. Ponto de Início (Círculo da mesma cor do sistema)
        push!(traces, scatter(
            x=[x_dados[1]], y=[y_dados[1]], 
            mode="markers", 
            name="Início ($nome_sistema)", 
            legendgroup="grupo_$idx_sol", 
            showlegend=false, # Esconde da legenda para não poluir
            marker=attr(size=8, color=cor_sistema, symbol="circle", line=attr(color="black", width=1))
        ))
        
        # C. Ponto Final (Quadrado da mesma cor do sistema)
        push!(traces, scatter(
            x=[x_dados[end]], y=[y_dados[end]], 
            mode="markers", 
            name="Fim ($nome_sistema)", 
            legendgroup="grupo_$idx_sol",
            showlegend=false,
            marker=attr(size=8, color=cor_sistema, symbol="square", line=attr(color="black", width=1))
        ))
    end
    
    # 3. O Alvo de Controle (Origem 0,0 - Adicionado apenas UMA vez)
    push!(traces, scatter(
        x=[0.0], y=[0.0], 
        mode="markers", 
        name="Origem", 
        marker=attr(size=12, color="black", symbol="x")
    ))
    
    # 4. Configuração Geométrica do Layout
    layout = Layout(
        title_text=titulo,
        title_x=0.5,
        xaxis_title="Estado x$(idx_x)",
        yaxis_title="Estado x$(idx_y)",
        width=700,
        height=700,
        showlegend=true,
        
        # Garante que seja um quadrado geométrico perfeito
        yaxis_scaleanchor="x",
        yaxis_scaleratio=1
    )
    
    # 5. Renderiza tudo
    fig = plot(traces, layout)
    display(fig)
end

function plotarRetratoFaseCompleto(dinamica, parametros=nothing;
                                   limiteEixo         = 3.0,
                                   densidadeSetas     = 18,
                                   raiosIniciais      = [0.5, 1.5, 2.5],
                                   trajetoriasPorAnel = 10,
                                   tempoMaximo        = 20.0,
                                   titulo             = "Retrato de Fase - Campo Vetorial")

    traces = GenericTrace[]
    escala = (2*limiteEixo / densidadeSetas) * 0.55

    # Adaptador: chama dinamica no formato DiffEq para obter (dx1, dx2)
    function avaliar(x1, x2)
        dx = zeros(2)
        dinamica(dx, [x1, x2], parametros, 0.0)
        return dx[1], dx[2]
    end

    # --- Campo Vetorial ---
    xs = range(-limiteEixo, limiteEixo, length=densidadeSetas)
    ys = range(-limiteEixo, limiteEixo, length=densidadeSetas)

    for x in xs, y in ys
        dx, dy = avaliar(x, y)
        mag = sqrt(dx^2 + dy^2)
        mag < 1e-8 && continue

        px = x + (dx/mag) * escala
        py = y + (dy/mag) * escala

        push!(traces, scatter(
            x=[x, px], y=[y, py],
            mode="lines",
            line=attr(width=1.2, color="rgba(100,130,220,0.5)"),
            showlegend=false, hoverinfo="skip"
        ))

        angulo_plotly = 90.0 - atan(dy, dx) * 180/π
        push!(traces, scatter(
            x=[px], y=[py],
            mode="markers",
            marker=attr(symbol="arrow", size=7, angle=angulo_plotly,
                        color="rgba(100,130,220,0.7)",
                        line=attr(width=1, color="rgba(100,130,220,0.7)")),
            showlegend=false, hoverinfo="skip"
        ))
    end

    push!(traces, scatter(x=[nothing], y=[nothing], mode="lines",
        line=attr(color="rgba(100,130,220,0.7)"), name="Campo vetorial"))

    # --- Trajetórias usando resolverSistema ---
    paleta = ["#2A5FCA","#1A8C60","#D05030","#8050C0","#1890A0","#C04060",
              "#5F8FD0","#40A860","#E07050","#A070D0","#40A0B8","#D06080"]

    idxCor = 1
    for raio in raiosIniciais
        angulos = range(0, 2π, length=trajetoriasPorAnel+1)[1:end-1]
        for a in angulos
            x0 = [raio * cos(a), raio * sin(a)]

            solucao = resolverSistema(dinamica, x0, (0.0, tempoMaximo), parametros)

            cor = paleta[mod1(idxCor, length(paleta))]
            push!(traces, scatter(
                x=[u[1] for u in solucao.u],
                y=[u[2] for u in solucao.u],
                mode="lines",
                line=attr(width=1.8, color=cor),
                legendgroup="traj",
                showlegend=(idxCor == 1),
                name="Trajetórias",
                text=["t = $(round(t, digits=2))" for t in solucao.t],
                hovertemplate="x₁ = %{x:.3f}<br>x₂ = %{y:.3f}<br>%{text}<extra></extra>"
            ))
            push!(traces, scatter(
                x=[solucao.u[1][1]], y=[solucao.u[1][2]],
                mode="markers",
                marker=attr(size=5, color=cor, symbol="circle"),
                legendgroup="traj",
                showlegend=false, hoverinfo="skip"
            ))
            idxCor += 1
        end
    end

    # Equilíbrio na origem
    push!(traces, scatter(
        x=[0.0], y=[0.0], mode="markers",
        marker=attr(size=10, color="red", symbol="x", line=attr(width=2, color="red")),
        name="Equilíbrio"
    ))

    layout = Layout(
        title_text=titulo, title_x=0.5,
        xaxis=attr(
            title="x₁",
            range=[-limiteEixo, limiteEixo],
            zeroline=true, zerolinecolor="#aaa", gridcolor="#ddd",
            tickmode="linear",
            dtick=round(limiteEixo/3, digits=1),
            showticklabels=true
        ),
        yaxis=attr(
            title="x₂",
            range=[-limiteEixo, limiteEixo],
            zeroline=true, zerolinecolor="#ddd", gridcolor="#ddd",
            tickmode="linear",
            dtick=round(limiteEixo/3, digits=1),
            showticklabels=true
        ),
        width=700, height=700, plot_bgcolor="#f8f7f4", showlegend=true
    )

    display(plot(traces, layout))
end
