"""
状态方程
"""
function pressure(;rho=0., e=0., gamma::Float64 = 1.4)
    return (gamma - 1.0) * rho * e
end

function pressure(w::Vector{Float64}, gamma::Float64)
    rho = w[1]
    e = w[end] / rho - 0.5 * MK.norm2(w[2:end-1] / rho)
    return (gamma - 1.0) * rho * e
end

function pressure_to_e(; rho::T = 0, p::T = 0, gamma::Float64 = 1.4) where T <: Real
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

function sound_speed(; rho::Float64 = 1, p::Float64 = 1, gamma::Float64 = 1.4)
    # if p < 0 || rho < 0
    #     # println(" p = ", p,"  rho = ",rho)
    #     error("sound_speed error")
    # end
    if rho == 0.
        return 0.
    elseif sign(rho) == sign(p)
        return  sqrt(gamma * p / rho)
    else
        return 0.
    end
end

function states2w(;rho::Float64=1.0,u::Array{Float64}=[0.], e::Float64=1.0)
    return [rho; rho * u; rho * (e+0.5*MK.norm2(u))]
end

function temperature(c::Cell, gamma::Float64)
    s = sound_speed(rho = c.rho, p = c.p, gamma = gamma)
    if isnan(s)
        println("-- temperature --")
        display(c)
        error("T error")
    end
    return s^2
end

function reynolds_number(rho::Float64, para::Dict)
    return rho * para["U0"] * para["L0"] / para["mu"]
end

function w2states!(c::Cell, gamma::Float64)
    if c.w[1] == 0.
        c.rho = 0.
        c.u = zeros(Float64, length(c.u))
        c.e = 0.
        c.p = 0. 
    else
        c.rho = c.w[1]
        c.u = c.w[2:end-1] / c.w[1]
        c.e = c.w[end] / c.w[1] - 0.5 * MK.norm2(c.u)
        c.p = pressure(rho = c.rho, e = c.e, gamma = gamma) 

        if c.rho < 0. && c.p > 0.
            c.rho = 0.
            c.u = zeros(Float64, length(c.u))
            c.e = 0.
            c.p = 0. 
            c.w = zeros(Float64, length(c.w))
        end
    end
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
    dim = length(W) - 2
    F = zeros(Float64, dim + 2)
    if W[1] != 0.
        rho, u, E = W[1], W[2:end-1]/W[1], W[end]/W[1]
        e = E - 0.5 * MK.norm2(u)
        p = pressure(rho = rho , e = e, gamma = gamma)
        # 数组的拼接非常耗时！！！改成下面这样可快一倍。
        F[1] = rho * u[1]
        F[2:end-1] = rho * u[axis] .* u
        F[1+axis] += p
        F[end] = (rho * E + p) * u[axis]
    end
    return F
end