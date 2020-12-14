
"""
MUSCL scheme
"""
function bi_interp!(axis::Int,WM2::Array{Float64,1},WM1::Array{Float64,1},WP1::Array{Float64,1},WP2::Array{Float64,1}, gamma::Float64)
    FM2,FM1,FP1,FP2 = w2f(axis,WM2,WM1,WP1,WP2, gamma::Float64)
    FL, FR = interp(FM2,FM1,FP1,FP2)
    WL, WR = interp(WM2,WM1,WP1,WP2)
    return FL,FR,WL,WR
end

function bi_interp!(axis::Int,WM1::Array{Float64,1},WP1::Array{Float64,1}, gamma::Float64)
    FM1, FP1 = w2f(axis,WM1,WP1, gamma::Float64)
    return FM1, FP1, WM1, WP1
end

function get_ausm_flux!(fL::Vector{Float64}, fR::Vector{Float64}, wL::Vector{Float64}, wR::Vector{Float64}, gamma::Float64; axis::Int = 0)
    len = length(fL)
    
    rhoL, uL, EL, pL, aL = get_flux_vars(wL, gamma)
    rhoR, uR, ER, pR, aR = get_flux_vars(wR, gamma)

    vel_L, vel_R = uL[axis], uR[axis]    

    a_int = max(0.5 * (aL + aR), 1.0)

    Ma_L, Ma_R = vel_L / a_int, vel_R / a_int
    
    Ma_bar_sq = (vel_L^2 + vel_R^2) / (2.0*a_int^2)

    fa = 1.0

    MMa_p4_L = get_MMa_4(Ma_L, fa,  1.0)
    MMa_m4_R = get_MMa_4(Ma_R, fa, -1.0)

    rho_int = (rhoL + rhoR) * 0.5

    if rho_int == 0.0 || fa == 0.0
        Ma_int = 0.0
    else
        Ma_int = MMa_p4_L + MMa_m4_R - AUSM_Kp/fa*max(1.0 - AUSM_sigma*Ma_bar_sq, 0.0)*(pR - pL)/(rho_int*a_int^2)
    end

    PP_p5_L = get_PP_5(Ma_L, fa,  1.0)
    PP_m5_R = get_PP_5(Ma_R, fa, -1.0)

    p_int = PP_p5_L*pL + PP_m5_R*pR - AUSM_Ku*PP_p5_L*PP_m5_R*(rhoL + rhoR)*fa*a_int*(vel_R - vel_L)

    vel_int = 0.5*(vel_L + vel_R)

    dm_int = a_int * Ma_int * (Ma_int > 0 ? rhoL : rhoR)
    
    # 绝对不要使用数组拼接 
    # Zha-Bilgen splitting
    Psi_int = ones(Float64, len)
    if dm_int > 0.0
        Psi_int[2:end-1] = uL
        Psi_int[end] = EL
    else
        Psi_int[2:end-1] = uR
        Psi_int[end] = ER
    end

    P_int = zeros(Float64, len)
    P_int[end] = p_int*vel_int
    P_int[1 + axis] = p_int

    f = dm_int * Psi_int + P_int

    return f
end

get_dflux!(reconst_scheme::String, ws::Array{Float64, 2}, para::Dict, axis::Int) = get_flux!(reconst_scheme, ws[:,1:4], para, axis = axis) - get_flux!(reconst_scheme, ws[:,2:5], para, axis = axis)

"""
'axis' stands for the space dimension of the flux.
"""
function get_flux!(reconst_scheme::String, ws::Array{Float64,2}, para::Dict; axis::Int = -10, flux_scheme::String = "AUSM")
    if reconst_scheme == "MUSCL"
        fL, fR, wL, wR = get_muscl_stencil_interp!(ws, para, axis = axis)
    elseif reconst_scheme == "WENO"
        error("undef reconst scheme")
    else
        error("undef reconst scheme")
    end

    f = get_stencil_flux!(fL, fR, wL, wR, para, axis = axis, flux_scheme = flux_scheme)

    # if f[1] != 0.0
    #     println("f = ", f)
    # end

    return f
end

function get_flux_vars(w::Vector{Float64}, gamma::Float64)
    if w[1] == 0.
        rho, u, E, p, a = 0., zeros(Float64, 3), 0., 0., 0.
    else
        rho = w[1]
        u = w[2:end-1] ./ w[1]
        E = w[end] / w[1]
        e = E - 0.5 * MK.norm2(u)
        p = pressure(w[1], e, gamma)
        # a = p < 0 ? 0.0 : sound_speed(rho = w[1], p = p, gamma = gamma)
        a = sound_speed(w[1], p, gamma)
    end
    return rho, u, E, p, a
end

function get_lf_flux!(FL::Vector{Float64},FR::Vector{Float64},WL::Vector{Float64},WR::Vector{Float64}, gamma::Float64)
    if WL[1] == 0.
        uL, EL, pL, WL, aL = zeros(Float64, 3), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uL = WL[2:end-1] ./ WL[1]
        EL = WL[end] / WL[1]
        eL = EL - 0.5 * MK.norm2(uL)
        pL = pressure(WL[1], eL, gamma)
        # if pL < 0.
        #     aL = 0.
        # else
            aL = sound_speed(WL[1], pL, gamma)
        # end
    end
    if WR[1] == 0.
        uR, ER, pR, WR, aR = zeros(Float64, 3), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uR = WR[2:end-1] ./ WR[1]
        ER = WR[4]/WR[1]
        eR = ER - 0.5 * MK.norm2(uR)
        pR = pressure(WR[1], eR, gamma)
        # if pR < 0.
        #     aR = 0.
        # else
            aR = sound_speed(WR[1], pR, gamma)
        # end
    end
    α = max(norm(uL)+aL, norm(uR)+aR)
    F = 0.5 .* ( FL .+ FR .- α .* ( WR .- WL ) )
    return F
end

function get_MMa_4(Ma::Float64, fa::Float64, s::Float64)
    if abs(Ma) >= 1
        MMa_1 = 0.5 * (Ma+s*abs(Ma))
        MMa_4 = MMa_1
    else
        MMa_2_p =   s*0.25*(Ma + s*1.0)^2
        MMa_2_m = - s*0.25*(Ma - s*1.0)^2
        MMa_4 = MMa_2_p*(1.0 - s*2.0*MMa_2_m)
    end
    return MMa_4
end

"""
MUSCL scheme
"""
function get_muscl_stencil_interp!(W::Array{Float64}, para::Dict; axis::Int = -10)
    # inviscous part
    if size(W,2) == 2
        FL, FR, WL, WR = bi_interp!(axis, W[:,1], W[:,2], para["gamma"])
    elseif size(W,2) == 4
        FL, FR, WL, WR = bi_interp!(axis, W[:,1], W[:,2], W[:,3], W[:,4], para["gamma"])
    else
        error( "Wrong size of W")
    end

    return FL, FR, WL, WR 
end

function get_PP_5(Ma::Float64, fa::Float64, s::Float64)
    if abs(Ma) >= 1 
        MMa_1 = 0.5*(Ma+s*abs(Ma))
        PP_5 = MMa_1 / Ma
    else
        MMa_2_p =   s*0.25*(Ma+s*1.0)^2
        MMa_2_m = - s*0.25*(Ma-s*1.0)^2
        alpha = 3.0/16.0 * (-4.0+5.0*fa^2)
        PP_5 = MMa_2_p * ((s*2.0 - Ma) - s*16.0*alpha*Ma*MMa_2_m)
    end
    return PP_5
end

function get_stencil_flux!(FL::Vector{Float64}, FR::Vector{Float64}, WL::Vector{Float64}, WR::Vector{Float64}, para::Dict; axis::Int = -10, flux_scheme::String = "AUSM")
    if flux_scheme == "LF"
        f = get_lf_flux!(FL, FR, WL, WR, para["gamma"])
    elseif flux_scheme == "AUSM"
        f = get_ausm_flux!(FL, FR, WL, WR, para["gamma"], axis = axis)
    else
        error("undef flux scheme")
    end
    return f 
end

"""
MUSCL scheme
"""
function interp(FM2::Array{Float64,1},FM1::Array{Float64,1},FP1::Array{Float64,1},FP2::Array{Float64,1})
    # MUSCL sample points: (FM2 - FM1 - pipe - FP1 - FP2 )
    n = length(FM2)
    FL = zeros(Float64,n)
    FR = zeros(Float64,n)
    for i = 1:n
        FL[i] = FM1[i] + 0.5 * minmod(FM1[i] - FM2[i], FP1[i] - FM1[i])
        FR[i] = FP1[i] - 0.5 * minmod(FP1[i] - FM1[i], FP2[i] - FP1[i])
    end
    return FL,FR
end

"""
minmod limiter
"""
function minmod(a::Float64, b::Float64)
    if abs(a) < abs(b)
        c = a
    elseif abs(a) > abs(b)
        c = b
    else
        c = 0.
    end
    return c
end