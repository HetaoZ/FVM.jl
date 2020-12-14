function pressure(rho, e, gamma)
    return (gamma - 1.0) * rho * e
end

function pressure(w::Vector{Float64}, gamma::Float64)
    return (gamma - 1.0) * w[1] * w[end] / w[1] - 0.5 * MK.norm2(w[2:end-1] / w[1])
end

function pressure_to_e(rho, p, gamma)
    if rho == 0.
        return 0.
    elseif sign(rho) == sign(p)
        return p / (gamma - 1.0) / rho 
    elseif rho > 0. && p <= 0.
        return p / (gamma - 1.0) / rho 
    else
        println("(rho, p) = ", (rho, p))
        error("wrong status")
    end
end

function sound_speed(rho, p, gamma)
    if rho == 0.
        return 0.
    elseif sign(rho) == sign(p)
        return  sqrt(gamma * p / rho)
    else
        return 0.
    end
end

function status_to_w(rho, u, e)
    return rho .* [1.0, u[1], u[2], u[3], (e+0.5*norm(u)^2)]
end

function temperature(rho, p, gamma)
    return sound_speed(rho, p, gamma)^2
end

function reynolds_number(rho::Float64, para::Dict)
    return rho * para["U0"] * para["L0"] / para["mu"]
end

function w_to_status(w, gamma::Float64)
    if w[1] == 0.
        rho = 0.
        u = zeros(Float64, 3)
        e = 0.
        p = 0. 
    else
        rho = w[1]
        u = w[2:4] / w[1]
        e = w[end] / w[1] - 0.5 * norm(u)^2
        p = pressure(rho, e, gamma) 

        if rho < 0. && p > 0.
            rho = 0.
            u = zeros(Float64, 3)
            e = 0.
            p = 0. 
        end
    end
    return rho, u, e, p
end

function w2f(axis::Int,WM2::Vector{Float64},WM1::Vector{Float64},WP1::Vector{Float64},WP2::Vector{Float64}, gamma::Float64)
    return w2f(axis, WM2, gamma::Float64), w2f(axis, WM1, gamma::Float64), w2f(axis, WP1, gamma::Float64), w2f(axis, WP2, gamma::Float64)
end

function w2f(axis::Int,WM1::Vector{Float64},WP1::Vector{Float64}, gamma::Float64)
    return w2f(axis, WM1, gamma::Float64), w2f(axis, WP1, gamma::Float64)
end

"""
Transform W to Flux in both dimensions
"""
function w2f(axis::Int, W::Vector{Float64}, gamma::Float64)
    F = zeros(Float64, 5)
    if W[1] != 0.
        rho, u, E = W[1], W[2:end-1]/W[1], W[end]/W[1]
        e = E - 0.5 * norm(u)^2
        p = pressure(rho, e, gamma)
        
        F[1] = rho * u[1]
        F[2:end-1] = rho * u[axis] .* u
        F[1+axis] += p
        F[end] = (rho * E + p) * u[axis]
    end
    return F
end