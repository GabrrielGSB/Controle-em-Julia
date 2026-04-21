# Controle/Exemplos/PID_no_pendulo_simples.jl

include("../SistemasLTI/Pêndulo.jl")
include("../controladores/PID.jl")
include("../Abstrações/MalhaFechada.jl")
include("../Ferramentas/ResoluçãoEDO.jl")
include("../Ferramentas/Visualização.jl")

# =========================================================================
# 1. Planta
# =========================================================================

pendulo = Pendulo(comprimento     = 1.0,
                  massa           = 10.0,
                  coefAtrito      = 0.1,
                  estadosIniciais = [1.0, 0.0])

# =========================================================================
# 2. Controlador
# =========================================================================

pid = PID(Kp = 200.0,
          Ki = 20.0,
          Kd = 10.0)

# =========================================================================
# 3. Malha Fechada
# =========================================================================

sys = conectar(pendulo, pid, π)      

# =========================================================================
# 4. Simulação
# =========================================================================

x0      = condicoesIniciais(sys, [1.0, 0.0])   
solucao = resolverSistema(sys, x0, (0.0, 10.0))

# =========================================================================
# 5. Visualização
# =========================================================================

plotarNoTempo(solucao, titulo="PID no Pêndulo Simples", estados=1:2)
gerarAnimacao(solucao, pendulo, 30)