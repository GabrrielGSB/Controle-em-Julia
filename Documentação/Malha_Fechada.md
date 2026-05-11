# Arquitetura e Implementação da Malha Fechada
Este documento detalha o funcionamento interno do módulo de Malha Fechada do FrameWork Controle. A arquitetura foi desenhada para ser modular, de alto desempenho ("Zero-Cost Abstraction") e perfeitamente integrada aos solvers numéricos do ecossistema SciML (como o OrdinaryDiffEq).

## 1. Visão Geral da Arquitetura
O framework adota uma abordagem moderna baseada em abstrações e Despacho Dinâmico, dividindo a responsabilidade do sistema em blocos isolados. O objetivo principal da estrutura de Malha Fechada é atuar como o mediador que faz a Planta (física) e o Controlador (lógica) conversarem no tempo certo a cada passo da integração numérica, sem que um precise conhecer os detalhes internos do outro. Essa separação resolve gargalos comuns em simulações de controle:

* **Fim do Acoplamento:** O controlador não precisa saber a quantidade total de estados da planta, operando apenas sobre a saída medida (o $y$).
  
* **Espaço de Estados Aumentado:** Estados internos do controlador (como a ação integral do PID) são tratados nativamente como novos estados do sistema diferencial acoplado.
  
* **Desempenho Otimizado:** Devido a ausência de variáveis globais e a especialização de tipos pelo compilador Julia é esperada velocidades de simulação comparáveis a C ou Fortran.

## 2. A Estrutura de Dados (struct MalhaFechada)
A conexão entre os elementos é materializada na struct ``MalhaFechada``. **Ela não deve ser instanciada diretamente** pelo usuário, mas sim através da função construtora ``conectar()``. 

A estrutura armazena:

* **Planta:** O objeto contendo os parâmetros físicos e a função de dinâmica.

* **Controlador:** O objeto contendo os ganhos e a lógica de controle.
  
* **Referencia:** O valor alvo desejado (setpoint) em formato de vetor.
  
* **``dim_planta`` e ``dim_ctrl``:** As dimensões dos estados físicos e de controle.
  
* **``idx_saida``:** Um vetor de inteiros (``Vector{Int}``) que define quais variáveis de estado da planta atuarão como as "saídas medidas" enviadas ao controlador.
  
### A Função conectar
A função conectar(; planta, controlador, referencia, idx_saida=1) atua como um empacotador. Ela inicializa a MalhaFechada e garante, por exemplo, que a referência e os índices de saída sejam formatados corretamente como vetores.

## 3. O Núcleo de Execução: O Functor
O grande "pulo do gato" da nossa arquitetura é a declaração do objeto ``MalhaFechada`` como um Functor (um objeto chamável).
Os solvers de equações diferenciais exigem que a função da dinâmica possua a assinatura estrita ``f(dx, x, p, t)``. Ao declarar a função diretamente no tipo ``MalhaFechada``, o objeto carrega consigo todos os parâmetros necessários (ganhos, massas, referências), evitando variáveis globais e desempacotamentos complexos.
O Functor é chamado milhares de vezes por segundo durante a simulação e executa os seguintes passos sequenciais:

### Passo 3.1: Fatiamento do Estado Aumentado (Sem Alocação)
O vetor de estado total ``x`` é dividido em duas partes usando macros ``@view`` para evitar alocação de memória:
* **``x_planta``**: Representa os estados físicos reais.
* **``x_ctrl``**: Representa a "memória" do controlador (ex: acúmulo da integral).
  
### Passo 3.2: Extração da Saída Medida ($y$)
O framework mapeia a equação de saída ($y = Cx$). O Functor extrai apenas as variáveis de interesse utilizando o campo de configuração:

```
y = x_planta[malha.idx_saida]
```
Isso encapsula a complexidade, entregando ao controlador apenas o que os sensores mediriam.

### Passo 3.3: Cálculo do Esforço de Controle ($u$)
O Functor aciona a interface do controlador, enviando apenas a saída filtrada ($y$), os estados internos do controlador, a referência e o tempo atual:Juliau = calcularSaida(malha.controlador, y, x_ctrl, malha.referencia, t)
Neste momento, a lógica matemática (como o cálculo de erro e multiplicação por ganhos PID) é resolvida.Passo 3.4: Evolução da Dinâmica da PlantaO esforço $u$ é injetado na física do sistema. A função de dinâmica específica da planta (ex: pendulo_invertido!) calcula as derivadas dos estados físicos (dx_planta).Passo 3.5: Evolução dos Estados Internos do ControladorSe o controlador possuir dinâmica interna (dim_ctrl > 0), o Functor aciona a atualização da memória do controlador:JuliaevoluirEstado!(malha.controlador, dx_ctrl, y, x_ctrl, malha.referencia, t)
Para um PID, por exemplo, é aqui que o integrador acumula o erro, definindo sua derivada (dx_ctrl) com base na distância entre a referência e a saída atual.4. O Fluxo de Vida na SimulaçãoMontagem: O usuário utiliza a função conectar() para unir Planta e Controlador, definindo os idx_saida.Condições Iniciais: A função condicoesIniciais() monta o vetor aumentado, onde os estados da planta recebem valores iniciais e a memória do controlador é inicializada em zero.Configuração do Solver: O objeto MalhaFechada é passado para o resolverSistema (que encapsula o ODEProblem). O pacote OrdinaryDiffEq recebe esse objeto e o enxerga como a própria equação diferencial.Integração: Durante o solve(), para cada passo de tempo $dt$, o solver invoca o Functor para descobrir a direção (dx) do sistema.Extração de Controle: Através de um SavingCallback configurado dentro do resolverSistema, o framework "espia" o cálculo de $u(t)$ a cada passo e o armazena em vetores históricos sem poluir a integração principal.