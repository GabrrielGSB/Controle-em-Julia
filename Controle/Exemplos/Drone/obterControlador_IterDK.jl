include("../../SistemasLTI/Drone.jl")
include("../../Não Linear/algoritmos/IteraçãoDK.jl")

f, g, vars = obterEqPolyDrone()

# =======================================================================================
#  OBTENÇÃO DO CONTROLADOR ALTITUDE
#  - Controla: Z (Altitude) e Vz (Velocidade no eixo Z)
#  - Entrada : U1 
    indices_descartados_altitude = [1, 2, 4, 5, 7, 8, 9, 10, 11, 12]
    subs_altitude = [vars[i] => 0.0 for i in indices_descartados_altitude]

    f_altitude    = [subs(f[3], subs_altitude...), subs(f[6], subs_altitude...)]
    g_altitude    = [g[3, 1], g[6, 1]]
    vars_altitude = [vars[3], vars[6]]

    Z, Vz = vars[3], vars[6]

    V0_altitude = 10.0*Z^2 + 2.0*Z*Vz + 5.0*Vz^2

    iteracaoDK(f_altitude, g_altitude, vars_altitude, V0_altitude,
               max_iter=1, nome_subsistema="Controlador de Altitude (Z)")
# =======================================================================================

# =======================================================================================
#  OBTENÇÃO DO CONTROLADOR DE ROLL (ROLAGEM)
#  Controla: Φ (Roll) e VΦ (Velocidade de Roll)
#  Entrada : U2 
    indices_descartados_roll = [1, 2, 3, 4, 5, 6, 8, 9, 11, 12]
    subs_roll = [vars[i] => 0.0 for i in indices_descartados_roll]

    f_roll    = [subs(f[7], subs_roll...), subs(f[10], subs_roll...)]
    g_roll    = [g[7, 2], g[10, 2]]
    vars_roll = [vars[7], vars[10]]

    Φ, VΦ = vars_roll[1], vars_roll[2]

    V0_roll = 10.0*Φ^2 + 2.0*Φ*VΦ + 5.0*VΦ^2

    iteracaoDK(f_roll, g_roll, vars_roll, V0_roll,
               max_iter=1, nome_subsistema="Controlador de Roll (Φ)")
# =======================================================================================

# =======================================================================================
#  OBTENÇÃO DO CONTROLADOR DE PITCH (ARFAGEM)
#  - Controla: θ (Pitch) e Vθ (Velocidade de Pitch)
#  - Entrada : U3 
    indices_descartados_pitch = [1, 2, 3, 4, 5, 6, 7, 9, 10, 12]
    subs_pitch = [vars[i] => 0.0 for i in indices_descartados_pitch]

    f_pitch    = [subs(f[8], subs_pitch...), subs(f[11], subs_pitch...)]
    g_pitch    = [g[8, 3], g[11, 3]]
    vars_pitch = [vars[8], vars[11]]

    θ, Vθ = vars_pitch[1], vars_pitch[2]

    V0_pitch = 10.0*θ^2 + 2.0*θ*Vθ + 5.0*Vθ^2

    iteracaoDK(f_pitch, g_pitch, vars_pitch, V0_pitch,
               max_iter=1, nome_subsistema="Controlador de Pitch (θ)")
# =======================================================================================

# =======================================================================================
#  OBTENÇÃO DO CONTROLADOR DE YAW (GUINADA)
#  - Controla: Ψ (Yaw) e VΨ (Velocidade de Yaw)
#  - Entrada : U4 
    indices_descartados_yaw = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11]
    subs_yaw = [vars[i] => 0.0 for i in indices_descartados_yaw]

    f_yaw    = [subs(f[9], subs_yaw...), subs(f[12], subs_yaw...)]
    g_yaw    = [g[9, 4], g[12, 4]]
    vars_yaw = [vars[9], vars[12]]

    Ψ, VΨ = vars_yaw[1], vars_yaw[2]

    V0_yaw = 10.0*Ψ^2 + 2.0*Ψ*VΨ + 5.0*VΨ^2

    iteracaoDK(f_yaw, g_yaw, vars_yaw, V0_yaw,
               max_iter=1, nome_subsistema="Controlador de Yaw (Ψ)")
# =======================================================================================