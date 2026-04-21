using LinearAlgebra

# 1. Estrutura de Parâmetros (agora recebendo U como entrada de controle virtual)
struct DroneParams
    massa           ::Float64
    gravidade       ::Float64
    Ixx             ::Float64
    Iyy             ::Float64
    Izz             ::Float64
    Ct              ::Float64 
    Cl              ::Float64 
    L               ::Float64 
    controleU       ::Vector{Float64} # Entradas virtuais [U1, U2, U3, U4]
    estadosIniciais ::NamedTuple
    variaveisEstado ::Tuple
    nomeSistema     ::String

    function DroneParams(; massa, gravidade=9.81, Ixx, Iyy, Izz, Ct, Cl, L, controleU=[0.0, 0.0, 0.0, 0.0], estadosIniciais)
        new(massa, gravidade, Ixx, Iyy, Izz, Ct, Cl, L, controleU, estadosIniciais, 
            ("X", "Y", "Z", "u", "v", "w", "ϕ", "θ", "ψ", "p", "q", "r"), 
            "Drone - Sistema Afim")
    end
end

# 2. Dinâmica Livre: f(x)
# Retorna um vetor 12x1 com as dinâmicas naturais (sem atuação dos motores)
function calcular_f(x, params)
    m, g = params.massa, params.gravidade
    Ixx, Iyy, Izz = params.Ixx, params.Iyy, params.Izz

    u, v, w_lin     = x[4], x[5], x[6]
    phi, theta, psi = x[7], x[8], x[9]
    p, q, r         = x[10], x[11], x[12]

    s_phi, c_phi = sin(phi), cos(phi)
    s_theta, c_theta, t_theta = sin(theta), cos(theta), tan(theta)
    s_psi, c_psi = sin(psi), cos(psi)

    fx = zeros(12)
    
    fx[1] = w_lin * (s_phi * s_psi + c_phi * c_psi * s_theta) - v * (c_phi * s_psi - c_psi * s_phi * s_theta) + u * c_psi * c_theta
    fx[2] = v * (c_phi * c_psi + s_phi * s_psi * s_theta) - w_lin * (c_psi * s_phi - c_phi * s_psi * s_theta) + u * c_theta * s_psi
    fx[3] = w_lin * c_phi * c_theta - u * s_theta + v * c_theta * s_phi
    
    fx[4] = r * v - q * w_lin + g * s_theta
    fx[5] = p * w_lin - r * u - g * c_theta * s_phi
    fx[6] = q * u - p * v - g * c_phi * c_theta
    
    fx[7] = p + q * s_phi * t_theta + r * c_phi * t_theta
    fx[8] = q * c_phi - r * s_phi
    fx[9] = r * (c_phi / c_theta) + q * (s_phi / c_theta)
    
    fx[10] = ((Iyy - Izz) / Ixx) * q * r
    fx[11] = ((Izz - Ixx) / Iyy) * p * r
    fx[12] = ((Ixx - Iyy) / Izz) * p * q

    return fx
end

# 3. Matriz de Atuação: g(x)
# Retorna uma matriz 12x4 que mapeia como os controles afetam os estados
function calcular_g(x, params)
    gx = zeros(12, 4)
    
    # Apenas os estados de velocidade Z e acelerações angulares recebem entrada
    gx[6, 1]  = 1.0 / params.massa
    gx[10, 2] = 1.0 / params.Ixx
    gx[11, 3] = 1.0 / params.Iyy
    gx[12, 4] = 1.0 / params.Izz
    
    return gx
end

# 4. Equação Diferencial Principal (Para o Solucionador ODE)
function drone_afim!(dx, x, params, t)
    # No futuro, o seu controle U virá de uma função, ex: U = lei_de_controle(x, params)
    # Por enquanto, pegamos o valor estático do parâmetro
    U = params.controleU 

    # dx = f(x) + g(x) * U
    fx = calcular_f(x, params)
    gx = calcular_g(x, params)
    
    # O operador .= atualiza o vetor dx inplace (muito mais rápido no Julia)
    dx .= fx + gx * U 
end