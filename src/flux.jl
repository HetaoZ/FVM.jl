
"""
MUSCL scheme
"""
function bi_interp!(axis::Int,WM2::Array{Float64,1},WM1::Array{Float64,1},WP1::Array{Float64,1},WP2::Array{Float64,1}, gamma::Float64)
    WL, WR = interp(WM2,WM1,WP1,WP2)
    
    # check if p > 0
    rhoL, uL, eL, pL = w_to_status(WL, gamma)
    rhoR, uR, eR, pR = w_to_status(WR, gamma)
    if pL <= 0
        WL = WM1
    end
    if pR <= 0
        WR = WP1
    end
    
    FL, FR = w2f(axis, WL, gamma), w2f(axis, WR, gamma)
    return FL, FR, WL, WR
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

    if a_int < 1e-14 || isnan(a_int)
        a_int = 1e14
    end

    Ma_L, Ma_R = vel_L / a_int, vel_R / a_int
    
    Ma_bar_sq = (vel_L^2 + vel_R^2) / (2.0*a_int^2)

    fa = 1.0

    MMa_p4_L = get_MMa_4(Ma_L, fa,  1.0)
    MMa_m4_R = get_MMa_4(Ma_R, fa, -1.0)

    rho_int = (rhoL + rhoR) * 0.5

    if rho_int < 1e-14 || fa < 1e-14
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

get_dflux!(id, reconst_scheme::String, ws::Array{Float64, 2}, para::Dict, axis::Int) = get_flux!(id, reconst_scheme, ws[:,1:4], para, axis = axis) - get_flux!(id, reconst_scheme, ws[:,2:5], para, axis = axis)

"""
'axis' stands for the space dimension of the flux.
"""
function get_flux!(id, reconst_scheme::String, ws::Array{Float64,2}, para::Dict; axis::Int = -10, flux_scheme::String = "AUSM")

    if reconst_scheme == "MUSCL"
        fL, fR, wL, wR = get_muscl_stencil_interp!(ws, para, axis = axis)
    elseif reconst_scheme == "WENO"
        error("undef reconst scheme")
    else
        error("undef reconst scheme")
    end

    wL = correct_cell_w(wL, para["gamma"], para["rho0"], para["u0"], para["e0"])
    wR = correct_cell_w(wR, para["gamma"], para["rho0"], para["u0"], para["e0"])

    rhoL, uL, eL, pL = FVM.w_to_status(wL, para["gamma"])
    rhoR, uR, eR, pR = FVM.w_to_status(wR, para["gamma"])

    if pL < 0
        wL = status_to_w(rhoL, uL,  para["e0"] * 1e-14)
    end
    if pR < 0
        wR = status_to_w(rhoR, uR,  para["e0"] * 1e-14)
    end
    
    # try
        f = get_stencil_flux!(fL, fR, wL, wR, para, axis = axis, flux_scheme = flux_scheme)

        return f
    # catch
    #     println("ws = ")
    #     display(ws)
    #     println()
    #     println(axis)
    #     println("wL = ",wL)
    #     println("wR = ",wR)
    #     println("fL = ",fL)
    #     println("fR = ",fR)
    #     error("error in w")
        
    # end
end

function get_flux_vars(w::Vector{Float64}, gamma::Float64)
    if w[1] < 1e-14
        rho, u, E, p, a = 0., zeros(Float64, 3), 0., 0., 0.
    else
        rho = w[1]
        u = w[2:end-1] ./ rho
        E = w[end] / rho
        e = E - 0.5 * norm(u)^2
        p = pressure(rho, e, gamma)
        if p < 0
        #     println("-- get_flux_vars: 1 --")
        #     println("p = ",p)
        #     println("w = ",w)
        #     error("")
            a = sound_speed(rho, 1e-14, gamma)
        else
            # a = p < 0 ? 0.0 : sound_speed(rho = w[1], p = p, gamma = gamma)
            a = sound_speed(rho, p, gamma)
        end
    end
    return rho, u, E, p, a
end

function get_lf_flux!(FL::Vector{Float64},FR::Vector{Float64},WL::Vector{Float64},WR::Vector{Float64}, gamma::Float64)
    if WL[1] < 1e-14
        uL, EL, pL, WL, aL = zeros(Float64, 3), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uL = WL[2:end-1] ./ WL[1]
        EL = WL[end] / WL[1]
        eL = EL - 0.5 * norm(uL)^2
        pL = pressure(WL[1], eL, gamma)
        # if pL < 0.
        #     aL = 0.
        # else
            aL = sound_speed(WL[1], pL, gamma)
        # end
    end
    if WR[1] < 1e-14
        uR, ER, pR, WR, aR = zeros(Float64, 3), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uR = WR[2:end-1] ./ WR[1]
        ER = WR[4]/WR[1]
        eR = ER - 0.5 * norm(uR)^2
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

# """
# rho-dominating MUSCL scheme
# """
# function interp(FM2::Array{Float64,1},FM1::Array{Float64,1},FP1::Array{Float64,1},FP2::Array{Float64,1})
#     # MUSCL sample points: (FM2 - FM1 - pipe - FP1 - FP2 )

#     if abs(FM1[1]-FM2[1]) > abs(FP1[1]-FM1[1])
#         selection = [0, 1]
#     else
#         selection = [1, 0]
#     end
#     FL = FM1 + 0.5 * selection' * [FM1-FM2, FP1-FM1]
    
#     if abs(FP1[1]-FM1[1]) > abs(FP2[1]-FP1[1])
#         selection = [0, 1]
#     else
#         selection = [1, 0]
#     end
#     FR = FP1 - 0.5 * selection' * [FP1-FM1, FP2-FP1]

    
#     return FL, FR
# end

"""
MUSCL scheme
"""
function interp(WM2::Array{Float64,1},WM1::Array{Float64,1},WP1::Array{Float64,1},WP2::Array{Float64,1})
    # MUSCL sample points: (WM2 - WM1 - pipe - WP1 - WP2 )
    
    # minmod limiter: old expression
    # FL =@. WM1 + 0.5 * minmod(WM1 - WM2, WP1 - WM1)
    # FR =@. WP1 - 0.5 * minmod(WP1 - WM1, WP2 - WP1)

    # mimod limiter: new expression
    rL = @. map(limited_r, WP1 - WM1, WM1 - WM2)
    rR = @. map(limited_r, WP2 - WP1, WP1 - WM1)

    opt = ("superbee", "van Leer mean", "van Leer", "van Albaba", "minmod")

    # compare dt at Step 20 and LF flux
    # LF:
    # (5,5) < (4,5) < (3,5) < (2,5) < (1,5) = 7.1e-4, small oscillations but unsymmetric
    # (5,4) = 6.7e-4, small oscillations and symmetric
    # (4,4) = 9.0e-5
    # 

    # (2,3) : no
    # (2,2): no
    #（3,4）:no

    opt1 = opt[5]
    opt2 = opt[5]

    WL = @. WM1 + 0.5 * (WM1 - WM2) * limiter(rL, opt1)
    WR = @. WP1 - 0.5 * (WP1 - WM1) * limiter(rR, opt1)

    WL[1] = WM1[1] + 0.5 * (WM1[1] - WM2[1]) * limiter(rL[1], opt2)
    WR[1] = WP1[1] - 0.5 * (WP1[1] - WM1[1]) * limiter(rR[1], opt2)

    

    return WL, WR
end

function limited_r(d1, d2)
    if d2 == 0.
        return 0.
    else
        return d1/ d2
    end
end

# -----------------------------------------------
# Limiters in high-to-low order of resolution of shock saves
# -----------------------------------------------
function limiter(r::Float64, name::String)
    if name == "superbee"
        return superbee_limiter(r)
    elseif name == "van Leer mean"
        return van_leer_mean_limiter(r)
    elseif name == "van Leer"
        return van_leer_limiter(r)
    elseif name == "van Albaba"
        return van_albaba_limiter(r)
    elseif name == "minmod"
        return minmod_limiter(r)
    else
        error("undef")
    end
end

"""
superbee limiter
"""
function superbee_limiter(r::Float64)
    return max(min(2*r,1), min(r,2))
end

"""
van Leer mean (double minmod) limiter
"""
function van_leer_mean_limiter(r::Float64)
    if r <= 0.
        return 0.
    else
        return minimum([2*r, 2, (1+r)/2])
    end
end

"""
van Leer limiter
"""
function van_leer_limiter(r::Float64)
    if 1 + r == 0.
        return 0.
    else
        return (r + abs(r)) / (1 + r)
    end
end

"""
van Albaba limiter
"""
function van_albaba_limiter(r::Float64)
    return (r^2+r)/(1+r^2)
end

"""
minmod limiter
"""
function minmod_limiter(r::Float64)
    if r <= 0.
        return 0.
    else
        return min(r, 1.)
    end
end

"""
minmod function
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