# =========================================================
# INCLUDES
    include("../../SistemasLTI/Drone.jl")
    include("../../Controladores/PID.jl")
    include("../../Abstrações/MalhaFechada.jl")
    include("../../Ferramentas/ResoluçãoEDO.jl")
    include("../../Ferramentas/Visualização.jl")
    include("../../Ferramentas/AnálisePerformance.jl")
# =========================================================

# =========================================================
# PARÂMETROS FÍSICOS
    m = 0.468
    g = 9.81

    drone = Drone(
        massa           = m,
        gravidade       = g,
        Ixx             = 4.856e-3,
        Iyy             = 4.856e-3,
        Izz             = 8.801e-3,
        Ct              = 2.980e-6,
        Cl              = 1.14e-7,
        L               = 0.225,
        estadosIniciais = zeros(12)
    )
# =========================================================

# =========================================================
# CONTROLADOR
    controlador = PID_MIMO(
        controladores = [PID(Kp = 2.0,  Ki = 0.5,  Kd = 1.0),   # Altitude (Z)
                         PID(Kp = 0.8,  Ki = 0.1,  Kd = 0.3),   # Roll     (Φ)
                         PID(Kp = 0.8,  Ki = 0.1,  Kd = 0.3),   # Pitch    (θ)
                         PID(Kp = 0.5,  Ki = 0.05, Kd = 0.2)],  # Yaw      (Ψ)
        
        ix_saida     = [3,  7,  8,  9],    # Z,  Φ,  θ,  Ψ
        idx_saida    = [6, 10, 11, 12],    # Vz, VΦ, Vθ, VΨ
        mapear_saida = (u, x) -> [
            clamp(m*g+u[1],  0.0,  2.0*m*g),  
            clamp(u[2],     -0.5,  0.5),       
            clamp(u[3],     -0.5,  0.5),       
            clamp(u[4],     -0.1,  0.1)        
        ]
    )
# =========================================================

# =========================================================
# SIMULAÇÃO
    ref = [1.0,   # Z_ref : Subir 1 metro
           0.0,   # Φ_ref : Roll  nivelado
           0.0,   # θ_ref : Pitch nivelado
           0.0]   # Ψ_ref : Yaw   nivelado

    sys = conectar(drone, controlador, ref)
    x0  = condicoesIniciais(sys, zeros(12))
    sol = resolverSistema(sys, x0, (0.0, 20.0),
                          resolucao       = 0.001,
                          salvar_controle = true)
# =========================================================

# =========================================================
# VISUALIZAÇÃO E ANÁLISE
    plotarNoTempo(sol,
                  titulo           = "Drone — Controle por PID_MIMO",
                  estados          = (3, 7, 8, 9),
                  mostrar_controle = true)

    # m_alt = analisarPerformance(sol.solucao,
    #                             referencia = ref[1],
    #                             idx_estado = 3)
    # imprimirRelatorio(m_alt, nome = "Canal Altitude (Z)")

    # m_roll = analisarPerformance(sol.solucao,
    #                              referencia = ref[2],
    #                              idx_estado = 7)
    # imprimirRelatorio(m_roll, nome = "Canal Roll (Φ)")

    # m_pitch = analisarPerformance(sol.solucao,
    #                               referencia = ref[3],
    #                               idx_estado = 8)
    # imprimirRelatorio(m_pitch, nome = "Canal Pitch (θ)")

    # m_yaw = analisarPerformance(sol.solucao,
    #                             referencia = ref[4],
    #                             idx_estado = 9)
    # imprimirRelatorio(m_yaw, nome = "Canal Yaw (Ψ)")
# =========================================================