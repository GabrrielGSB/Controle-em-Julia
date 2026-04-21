using DynamicPolynomials

"""
    Recebe os polinômios do sistema e retorna uma função rápida no formato `(dx, x, p, t)` 
    pronta para ser usada no ODEProblem.
"""
function extrairSistemaEDO(f_x, vars; g_x=nothing, u_x=nothing)
    if g_x !== nothing && u_x !== nothing
        dinamica_total = f_x + g_x .* Ref(u_x)
    else
        dinamica_total = f_x
    end

    return function sistema!(dx, x, p, t)
        for i in 1:length(dx)
            dx[i] = dinamica_total[i](vars => x)
        end
    end
end