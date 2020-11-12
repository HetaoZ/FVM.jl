
module FVM
# --------------------------------------
# 接口
using PyPlot # "应该把后处理分割出去"
using DelimitedFiles, Printf, Distributed, DistributedArrays, LinearAlgebra, Statistics
using MathKits
const MK = MathKits
export Fluid, copy!, GCell, fill_fluid!, set_boundaries!, after_shock
#--------------------------------------
# 常数 
const OUTPUTDATA = true
const AUSM_Kp = 0.25
const AUSM_sigma = 1
const AUSM_Ku = 0.75
const MAX_DIM = 3
const DEBUG_MODE = true
const NP = 1  # DEBUG_TEST
const TOL_STEP = 1.e-10
const FIGHEIGHT = 8
const FRAME_BASE = 1000000
const RK_COEFF = [1.0 0.75 1/3;
                  0.0 0.25 2/3;
                  1.0 0.25 2/3]
# --------------------------------------

# -----------------------------------------------------------------------------------
# 宏


# 函数

function fvm_advance!(f::Fluid, dt::Float64)
    
    # println("-- advance_fluid: 1 --")
    # showfield!(f.cells, "rho", 13:20)

    backup_w!(f)

    # println("-- advance_fluid: 2 --")
    # showfield!(f.cells, "rho", 13:20)

    # @warn("rk is limited")
    for rk = 1:3
        # println("-- rk = ",rk, " --")

        # println("-- advance_fluid: 3 --")
        # showfield!(f.cells, "flux", 13:20)

        update_rhs!(f)

        # println("-- advance_fluid: 4 --")
        # showfield!(f.cells, "rhs", 13:20)

        update_cells!(f, rk, dt)

        # println("-- advance_fluid: 5 --")
        # showfield!(f.cells, "rho", 13:20)

        # 下面这种写法是错的，因为wb不会随着rk更新。
        # update_fluxes_and_cells!(f, rk = rk, dt = dt)
        update_boundaries!(f)

        # println("-- advance_fluid: 6 --")
        # showfield!(f.cells, "rho", 13:20)
    end

end

# """
# Return the mass density, the pressure, and the velocity (relative to the ground) of the fluid after a shock wave.
# """
# function after_shock(; p::T = 0, rho::T = 0, Mach::T = 0, gamma::T = 1.4) where T <: Real
#     @assert Mach >= 1
#     cs2 = gamma/(gamma+1) * p / rho * (2 + (gamma-1)*Mach^2)
#     U1  = Mach * sqrt(gamma*p/rho)
#     U2  = cs2 / U1
#     p2  = p - rho * U1 * (U2-U1)
#     rho2  = 1 / (1/(rho*U1)*(U2+U1) - 1/rho)
#     return rho2, p2, -(U2 - U1)
# end

"""
向sdir方向运动的激波，波前是区域1，波后是区域2，激波本身记为s，激波马赫数Ms=激波的地面速度us的绝对值除以波前声速c1。

算例：Monasse(2012) Example 2

mshock(1.01e5,1.2,0,1.21,1.4,-1) = (1.630779226806516, 155686.45, -109.71829086299334)
"""
function after_shock(p1,ρ1,u1,Ms,γ,sdir)
    @assert Ms >= 1
    c1 = sound_speed(rho = ρ1, p = p1, gamma = γ)
    us = Ms * c1 * sdir
    U1 = u1 - us
    M1 = U1 / c1
    r = (γ+1)*M1^2 / (2 + (γ-1)*M1^2)
    U2 = U1 / r
    u2 = U2 + us
    ρ2 = ρ1 * r
    p2 = p1 * (1 + 2*γ/(γ+1) * (M1^2-1))
    # c2 = sp(γ,p2,ρ2)
    # M2 = abs(U2 / c2)
    # M2_predicted = sqrt( (1+(γ-1)/2*M1^2) / (γ*M1^2-(γ-1)/2) )
    # println(M2)
    # println(M2_predicted)
    return ρ2, p2, u2
end

function backup_w!(f::Fluid)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                c.wb = copy(c.w)
            end
        end
    end   
end

function testbackup_w!(f::Fluid)
    @sync for pid in workers()
        @spawnat pid begin
            @inbounds @simd for c in localpart(f.cells)
                c.wb = copy(c.w)
            end
        end
    end   
end

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

function check_cell_status!(c::GCell)
    nan = false
    if isnan(sum(c.w)) || isnan(sum(c.rhs)) || isnan(sum(c.wb))
        nan = true
    end

    if nan
        println("-- c is NaN --")
        display(c)
        error("NaN error")
    end
end

function check_all_cells_status!(cells::DArray)
    for c in cells
        check_cell_status!(c)
    end
end

function check_all_cells_status!(f::Fluid)
    for c in f.cells
        check_cell_status!(c)
    end
end

function check_all_cells_particles_status!(f::Fluid)
    check_all_cells_status!(f.cells)
    check_all_particles_status!(f.particles)
end

function check_all_particles_status!(particles::Array{Particle})
    for p in particles
        check_particle_status!(p)
    end
end

function check_all_particles_status!(f::Fluid)
    for p in f.particles
        check_particle_status!(p)
    end
end

function check_conservativity!(f::Fluid)
    @warn "This may take much time!"
    mass = 0.
    for c in f.cells
        if MK.betweeneq(c.x, f.point1, f.point2)
            mass += c.rho * prod(f.d) 
        end
    end
    if !f.total_is_summed
        f.total_mass = mass
        f.total_is_summed = true
    end
    println("-- [ mass = ", mass, " | δ = ", (mass-f.total_mass)/f.total_mass, " ] --")    
end

function check_mass!(f::Fluid; point1::Array{Float64} = f.point1, point2::Array{Float64} = f.point2)
    mass = 0.
    for c in f.cells
        if MK.betweeneq(c.x, point1, point2)
            mass += c.rho * prod(f.d) 
        end
    end
    if !f.total_is_summed
        f.total_mass = mass
        f.total_is_summed = true
    end
    
    return mass
end

function check_mass!(particles::Array{Particle})
    mass = 0.
    for p in particles
        mass += p.m
    end
    return mass
end

function check_particle_status!(p::Particle)
    nan = false
    if isnan(sum(p.m)) || isnan(sum(p.mu)) || isnan(sum(p.mE))
        nan = true
    end
    if isnan(sum(p.x)) || isnan(sum(p.dx)) || isnan(sum(p.force_to_boundary))
        nan = true
    end
    if isnan(sum(p.boundary_id)) || isnan(sum(p.boundary_ratio))
        nan = true
    end

    if nan
        println("-- p is NaN --")
        display(p)
        error("NaN error")
    end    
end

function Base.copy!(f::Fluid)
    f1 = deepcopy(f)
    f1.cells = deepcopy(f.cells)
    # f1 = Fluid(f.dim)
    # f1.point1 = copy(f.point1)
    # f1.point2 = copy(f.point2)
    # f1.nmesh = copy(f.nmesh)
    # f1.d = copy(f.d)
    # f1.dist = copy(f.dist)
    # f1.boundaries = copy(f.boundaries)
    # f1.ng = copy(f.ng)
    # f1.ndiv = copy(f.ndiv)
    # f1.cells = copy!(f.cells, f.ndiv, f.dist)
    # f1.particles = copy!(f.particles)
    # f1.constants = f.constants
    # f1.background_is_filled = f.background_is_filled
    # f1.total_mass = f.total_mass
    # f1.total_is_summed = f.total_is_summed
    # f1.consider_vis_item = f.consider_vis_item
    # f1.exclude_particles = f.exclude_particles
    # f1.reconst_scheme = f.reconst_scheme
    # f1.flux_scheme = f.flux_scheme
    return f1
end

function Base.copy!(particles::Array{Particle})
    ps = Array{Particle}(undef, size(particles))
    for k in eachindex(particles)
        ps[k] = copy!(particles[k])
    end
    return ps
end

function Base.copy!(particle::Particle)
    p = Particle(particle.dim)
    p.m = particle.m
    p.u = copy(particle.u)
    p.ek = particle.ek
    p.e = particle.e
    p.E = particle.E
    p.x = copy(particle.x)
    p.dx = copy(particle.dx)
    p.V  = particle.V
    p.r  = particle.r
    p.force_to_boundary = copy(particle.force_to_boundary)
    p.boundary_id = particle.boundary_id
    p.boundary_ratio = particle.boundary_ratio
    return p
end

function Base.copy!(c::GCell) 
    cell = GCell(length(c.u))
    # 必须复制坐标
    cell.rho = c.rho
    cell.u = copy(c.u)
    cell.e = c.e
    cell.p = c.p
    cell.w = copy(c.w)
    cell.wb = copy(c.wb)   
    cell.rhs = copy(c.rhs)
    cell.i = copy(c.i)
    cell.x = copy(c.x)
    return cell
end


function Base.copy!(a::Array{Vector{Union{Nothing, Float64}}})
    b = Array{Vector{Union{Nothing, Float64}}}(undef, size(a))
    for i in eachindex(a)
        b[i] = copy(a[i])
    end
    return b
end

function Base.copy!(a::Vector{Union{Nothing, Float64}})
    b = Vector{Union{Nothing, Float64}}(undef, size(a))
    for i in eachindex(a)
        b[i] = a[i]
    end
    return b
end

function Base.copy!(a::Array{GCell}) 
    b = Array{GCell}(undef, size(a))
    for k in eachindex(b)
        b[k] = copy!(a[k])
    end
    return b
end

"""
这个函数如果写的方式不对，就会相当耗时。
"""
function Base.copy!(a::DArray{GCell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
    dim = length(a[1].u)
    tmp = Array{GCell}(undef, size(a))
    for i in eachindex(tmp)
        tmp[i] = GCell(dim)
    end
    b = distribute(tmp, procs = workers(), dist = dist)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(a)
                copy_complete_cell_to!(c, b[Tuple(c.i)...])
            end
        end
    end
    return b
end

"""
这个函数如果写的方式不对，就会相当耗时。
"""
function testcopy!(a::DArray{GCell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
    dim = a[1].dim
    tmp = Array{GCell}(undef, size(a))
    for i in eachindex(tmp)
        tmp[i] = GCell(dim)
    end
    b = distribute(tmp, procs = workers(), dist = dist)
    if dim == 1
        @sync for pid in workers()
            @spawnat pid begin
                for c in localpart(a)
                    copy_complete_cell_to!(c, b[Tuple(c.i)...])
                end 
            end
        end
    elseif dim == 2
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(a)
                for i in inds[1]
                    for j in inds[2]
                        copy_complete_cell_to!(a[i, j], b[i, j])
                    end
                end 
            end
        end
    elseif dim == 3
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(a)
                for i in inds[1]
                    for j in inds[2]
                        for k in inds[3]
                            copy_complete_cell_to!(a[i, j, k], b[i, j, k])
                        end
                    end
                end 
            end
        end
    else
        error("undef dim")
    end
    return b
end

# function Base.copy!(a::DArray{GCell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
#     tmp = Array{GCell}(undef, size(a))
#     @time for i in eachindex(a)
#         c = a[i]
#         tmp[i] = copy!(c)
#     end
    
#     b = distribute(tmp, procs = workers(), dist = dist)
#     return b
# end

function Base.copy!(a::SubArray{GCell,N,P,I,L} where L where I where P<:DArray where N)
    b = Array{GCell}(undef, size(a))
    for k in eachindex(b)
        b[k] = copy!(a[k])
    end
    return b
end

function copy_to!(c::GCell, target::GCell)
    target.rho = c.rho
    target.u = copy(c.u)
    target.e = c.e
    target.p = c.p
    target.w = copy(c.w)
    target.wb = copy(c.wb)
end

function copy_to!(a::Array{GCell}, target::Array{GCell}) 
    @assert size(a) == size(target)
    for k in eachindex(target)
        copy_to!(a[k], target[k])
    end
end
function copy_to!(a::SubArray, target::SubArray) 
    @assert size(a) == size(target)
    for k in eachindex(target)
        copy_to!(a[k], target[k])
    end
end
function copy_to!(a::SubArray, target::Array{GCell}) 
    @assert size(a) == size(target)
    for k in CartesianIndices(target)
        # println("-- copy_to: 1 --")
        # println(k)
        # println(typeof(a))
        # println(typeof(target))
        # println(typeof(a[k]))
        # println(typeof(target[k]))
        copy_to!(a[k], target[k])
    end
end

function copy_complete_cell_to!(c::GCell, target::GCell)
    target.rho = c.rho
    target.u = copy(c.u)
    target.e = c.e
    target.p = c.p
    target.w = copy(c.w)
    target.wb = copy(c.wb)
    target.rhs = copy(c.rhs)
    target.i = copy(c.i)
    target.x = copy(c.x)
end

"""
在密度相对较低的单元，数值伪振荡可能造成负质量或负压力。对其进行保正性修正，会造成轻微不守恒，总质量会轻微增加。
"""
function correct_cell_w!(c::GCell, gamma)
    if sign.([c.w[1], c.w[end], pressure(c.w, gamma)]) == [-1.0, -1.0, -1.0]
        # 允许负质量单元存在。
    elseif c.w[1] < 0 || c.w[end] < 0 || pressure(c.w, gamma) < 0
        # println("-- correct_cell_w: 1 --")
        # println("c.i = ", c.i)

        c.w = zeros(Float64, length(c.w))
    end
end

function fill_fluid!(f::Fluid, cell::GCell)
    f.cells = distribute(reshape([copy!(cell) for i in 1:prod(f.nmesh .+ f.ng*2)], Tuple(f.nmesh .+ f.ng*2)...), procs = workers(), dist = f.dist)
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.cells)
            bias = [inds[k][1] - 1 for k = 1:f.dim]
            local_inds = Tuple([inds[k] .- bias[k] for k in 1:f.dim])
            local_cells = getindex(localpart(f.cells), local_inds...)
            for i in CartesianIndices(local_cells)
                c = local_cells[i]
                c.i = [i[k] + bias[k] for k = 1:f.dim]
                c.x = f.d .* (c.i .- f.ng .- 0.5) .+ f.point1
            end
        end
    end
    f.background_is_filled = true  
end

function fill_fluid!(f::Fluid, cell::GCell, point1::Array, point2::Array)
    if !f.background_is_filled
        error("Please fill the background before this calling.")
    else
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(f.cells)
                bias = [inds[k][1] - 1 for k = 1:f.dim]
                local_inds = Tuple([inds[k] .- bias[k] for k in 1:f.dim])
                local_cells = getindex(localpart(f.cells), local_inds...)
                for c in local_cells
                    if MK.between(c.x, point1, point2)
                        copy_to!(cell, c)
                    end    
                end
            end
        end        
    end
    
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

get_dflux!(reconst_scheme::String, ws::Array{Float64, 2}, constants::Dict, axis::Int) = get_flux!(reconst_scheme, ws[:,1:4], constants, axis = axis) - get_flux!(reconst_scheme, ws[:,2:5], constants, axis = axis)

"""
这里的axis表示通量所在维度。
"""
function get_flux!(reconst_scheme::String, ws::Array{Float64,2}, constants::Dict; axis::Int = -10, flux_scheme::String = "AUSM")
    if reconst_scheme == "MUSCL"
        fL, fR, wL, wR = get_muscl_stencil_interp!(ws, constants, axis = axis)
    elseif reconst_scheme == "WENO"
        error("undef reconst scheme")
    else
        error("undef reconst scheme")
    end

    f = get_stencil_flux!(fL, fR, wL, wR, constants, axis = axis, flux_scheme = flux_scheme)


    return f
end

function get_flux_vars(w::Vector{Float64}, gamma::Float64)
    if w[1] == 0.
        rho, u, E, p, a = 0., zeros(Float64, length(w) - 2), 0., 0., 0.
    else
        rho = w[1]
        u = w[2:end-1] ./ w[1]
        E = w[end] / w[1]
        e = E - 0.5 * MK.norm2(u)
        p = pressure(rho = w[1], e = e, gamma = gamma)
        # a = p < 0 ? 0.0 : sound_speed(rho = w[1], p = p, gamma = gamma)
        a = sound_speed(rho = w[1], p = p, gamma = gamma)
    end
    return rho, u, E, p, a
end

function get_lf_flux!(FL::Vector{Float64},FR::Vector{Float64},WL::Vector{Float64},WR::Vector{Float64}, gamma::Float64)
    if WL[1] == 0.
        uL, EL, pL, WL, aL = zeros(Float64, length(WL) - 2), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uL = WL[2:end-1] ./ WL[1]
        EL = WL[end] / WL[1]
        eL = EL - 0.5 * MK.norm2(uL)
        pL = pressure(rho = WL[1], e = eL, gamma = gamma)
        # if pL < 0.
        #     aL = 0.
        # else
            aL = sound_speed(rho = WL[1], p = pL, gamma = gamma)
        # end
    end
    if WR[1] == 0.
        uR, ER, pR, WR, aR = zeros(Float64, length(WL) - 2), 0., 0., zeros(Float64, length(WL)), 0.
    else
        uR = WR[2:end-1] ./ WR[1]
        ER = WR[4]/WR[1]
        eR = ER - 0.5 * MK.norm2(uR)
        pR = pressure(rho = WR[1], e = eR, gamma = gamma)
        # if pR < 0.
        #     aR = 0.
        # else
            aR = sound_speed(rho = WR[1], p = pR, gamma = gamma)
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
function get_muscl_stencil_interp!(W::Array{Float64}, constants::Dict; axis::Int = -10)
    # 非粘性部分
    if size(W,2) == 2
        FL, FR, WL, WR = bi_interp!(axis, W[:,1], W[:,2], constants["gamma"])
    elseif size(W,2) == 4
        FL, FR, WL, WR = bi_interp!(axis, W[:,1], W[:,2], W[:,3], W[:,4], constants["gamma"])
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

function get_stencil_flux!(FL::Vector{Float64}, FR::Vector{Float64}, WL::Vector{Float64}, WR::Vector{Float64}, constants::Dict; axis::Int = -10, flux_scheme::String = "AUSM")
    if flux_scheme == "LF"
        f = get_lf_flux!(FL, FR, WL, WR, constants["gamma"])
    elseif flux_scheme == "AUSM"
        f = get_ausm_flux!(FL, FR, WL, WR, constants["gamma"], axis = axis)
    else
        error("undef flux scheme")
    end
    return f 
end

get_sxx(hx::Float64, u::Vector{Float64}, hy::Float64, v::Vector{Float64}, mu::Float64) = 2/3*mu*(2*(u[3] - u[1])/(2*hx) - (v[3] - v[1])/(2*hy))

get_sxy(hy::Float64, u::Vector{Float64}, hx::Float64, v::Vector{Float64}, mu::Float64) = mu*((u[3] - u[1])/(2*hy) + (v[3] - v[1])/(2*hx))

get_syy(hx::Float64, u::Vector{Float64}, hy::Float64, v::Vector{Float64}, mu::Float64) = 2/3*mu*(- (u[3] - u[1])/(2*hx) + 2*(v[3] - v[1])/(2*hy))

function get_vis_item!(u::Array{Float64}, v::Array{Float64}, T::Array{Float64}, constants::Dict, d::Vector, rho::Float64)
    if rho == 0.0
        return zeros(Float64, length(d) + 2)
    end
    Re = reynolds_number(rho, constants)
    mu = constants["mu"]
    beta = - mu/(constants["Pr"]*(constants["gamma"]-1.0))
    if length(d) == 2
        hx, hy = d[1], d[2]

        # 现在在5x5的点阵里考虑. 编号1~5，中心是3.
        sxx_43 = get_sxx(hx, u[3:5, 3], hy, v[4, 2:4], mu)
        sxx_23 = get_sxx(hx, u[1:3, 3], hy, v[2, 2:4], mu)
    
        sxy_43 = get_sxy(hy, u[4, 2:4], hx, v[3:5, 3], mu) # 注意应力张量对称
        sxy_23 = get_sxy(hy, u[2, 2:4], hx, v[1:3, 3], mu) # 注意应力张量对称
    
        sxy_34 = get_sxy(hy, u[3, 3:5], hx, v[2:4, 4], mu)
        sxy_32 = get_sxy(hy, u[3, 1:3], hx, v[2:4, 2], mu)
    
        syy_34 = get_syy(hx, u[2:4, 4], hy, v[3, 3:5], mu)
        syy_32 = get_syy(hx, u[2:4, 2], hy, v[3, 1:3], mu)
    
        qx_43 = beta * (T[5, 3] - T[3, 3])/(2*hx)
        qx_23 = beta * (T[3, 3] - T[1, 3])/(2*hx)
        
        qy_34 = beta * (T[3, 5] - T[3, 3])/(2*hy)
        qy_32 = beta * (T[3, 3] - T[3, 1])/(2*hy)
    
        bx_43 = u[4, 3]*sxx_43 + v[4, 3]*sxy_43 - qx_43
        bx_23 = u[2, 3]*sxx_23 + v[2, 3]*sxy_23 - qx_23
    
        by_34 = u[3, 4]*sxy_34 + v[3, 4]*syy_34 - qy_34
        by_32 = u[3, 2]*sxy_32 + v[3, 2]*syy_32 - qy_32
        
        vis_item_x = [0.0, sxx_43 - sxx_23, sxy_43 - sxy_23, bx_43 - bx_23] ./ (2*hx) 
        vis_item_y = [0.0, sxy_34 - sxy_32, syy_34 - syy_32, by_34 - by_32] ./ (2*hy) 
    
        vis_item = (vis_item_x .+ vis_item_y) ./ Re
        if any(map(isnan, vis_item))
            println("-- get_vis_item: 0 --")
            display(vis_item_x)
            display(vis_item_y)
            display(Re)
            display([sxx_43 sxx_23 sxy_43 sxy_23 bx_43 bx_23 hx;
                     sxy_34 sxy_32 syy_34 syy_32 by_34 by_32 hy])
            error("get_vis_item error")
        end
    else
        error("get_vis_item error")
    end
    return vis_item
end

"""
gcells是要被修改的边界外单元。cells是作为参照的边界内单元。
"""
function ghost_boundary_cell!(cells::SubArray, gcells::Array{GCell}, boundary::String; axis::Int = -10)
    copy_to!(cells, gcells)
    if boundary == "refl"
        factor = - 1
    else
        factor = 1
    end

    for k in eachindex(cells)

        gcells[k].u[axis] *= factor
        gcells[k].w[axis + 1] *= factor
    end 
end

function index_to_pid(k::Union{Array{Int}, CartesianIndex}; nzone::Array{Int} = [1], dist::Array{Int} = [1])
    K = cld.(Tuple(k), nzone ./ dist)
    pid = workers()[Int(K[1] + [K[i] - 1 for i in 2:length(K)]' * [prod(dist[1:i]) for i in 1:length(K)-1])]
    return pid
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

function length_of_tuples(r::Tuple)
    s = 1
    for k in r
        s *= length(k)
    end
    return s
end

function meshxy(x::Array, y::Array; dim::Int = 1)
    m=maximum(size(x))
    n=maximum(size(y))
    X=Array{Float64}(undef,m,n)
    Y=Array{Float64}(undef,m,n)
    for i=1:m,j=1:n
        X[i,j]=x[i]
        Y[i,j]=y[j]
    end
    if dim == 1
        return reshape(X,m*n,1),reshape(Y,m*n,1)
    elseif dim == 2
        return X,Y
    else
        error( "The dimension must be <= 2")
    end
end

"""
minmod limiter
"""
function minmod(a::Float64, b::Float64)
    if abs(a) < abs(b)
        return a
    elseif abs(a) > abs(b)
        return b
    else
        return 0.
    end
end

function output!(frame::Int, time::Float64; filepath::String = "../outputdata/")
    if frame == 0
        open(filepath*"time.txt","w") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    else
        open(filepath*"time.txt","a") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    end
end

"""
axis表示输出变量的方向。
"""
function output!(f::Fluid; varname::String = "rho", axis::Int = 1, frame::Int = 1, filepath::String = "../outputdata/")
    if varname == "mesh"
        for k in 1:f.dim
            open(filepath*"fluid_x"*string(k)*".txt","w") do file
                writedlm(file, [f.d[k]*(i-0.5-f.ng)+f.point1[k]  for  i=1:f.nmesh[k]+f.ng*2])
            end
        end        

    else
        # essential vars
        for essential_varname in ("rho", "p")
            data = Array{typeof(getfield(f.cells[1],Symbol(essential_varname)))}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k],Symbol(essential_varname))
            end      
            open(filepath*"fluid_"*essential_varname*"_"*string(frame+FRAME_BASE)*".txt","w") do file
                writedlm(file, data)
            end           
        end           
        # reviewed vars
        if varname in ("e",)
            data = Array{typeof(getfield(f.cells[1],Symbol(varname)))}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k],Symbol(varname))
            end
        end
        if varname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end      
            open(filepath*"fluid_"*varname*"_"*string(frame+FRAME_BASE)*".txt","w") do file
                writedlm(file, data)
            end      
        end

    end
end

function output_plot_func(x::Array{Float64}, data::Array{Float64}; levels::Vector = [0, 1, 1])
    figure(1, figsize=(8, 6))
    c = plot(x, data)

    # ref
    A = readdlm("../../../CommonData/monasse2012-piston-p.txt")
    plot(A[:,1], A[:,2], color = "r")
    
    ylim(levels[1], levels[2])
end

function output_plot_func(x::Array{Float64}, y::Array{Float64}, data::Array{Float64}; levels::Vector = [0, 1, 1], cmap::String = "jet", cbar::Bool = true)
    scale = FIGHEIGHT / (y[end] - y[1])
    figure(1, figsize=((x[end] - x[1]) * scale + 4 * cbar, FIGHEIGHT))
    X, Y = meshxy(x, y, dim = 2)
    if levels[2] <= levels[1] || levels[3] < 2
        levels = (minimum(data), maximum(data), 30)
    end
    c = contourf(X, Y, data, levels=LinRange(Tuple(levels)...), cmap=cmap);
    c1 = contour(X, Y, data, levels=LinRange(Tuple(levels)...));
    # c2 = scatter(meshxy(x, y, dim = 1)..., marker = "+");
    xlim([x[3], x[end-2]]); ylim([y[3], y[end-2]])
    if cbar 
        colorbar(c)
    end
end

function outputfig(frame::Int; varname::String = "rho", filepath::String = "outputdata/", figpath::String = "outputfig/",  levels::Vector = [0, 1, 1], cmap::String = "jet", cbar::Bool = true, plotdim::String = "2D")
    if plotdim == "1D"
        x = read_mesh(1, frame, varname = varname, filepath = filepath)
        data = read_frame(frame, varname = varname, filepath = filepath)
        ioff()
        output_plot_func(x[1], data, levels=levels)
    elseif plotdim in ("2D", "2D-x", "2D-y")
        x, y = read_mesh(2, frame, varname = varname, filepath = filepath)
        data = read_frame(frame, varname = varname, filepath = filepath)
        ioff()
        if plotdim == "2D-x"
            output_plot_func(x, data[:, floor(Int, size(data, 2)/2)], levels=levels)
        elseif plotdim == "2D-y"
            output_plot_func(y, data[floor(Int, size(data, 1)/2), :], levels=levels)  
        else
            output_plot_func(x, y, data, levels = levels, cmap = cmap, cbar = cbar)    
        end 
    else
        error("undef dim")       
    end
    savefig(figpath*varname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end

function outputfig(frames::AbstractRange; varname::String = "rho", filepath::String = "outputdata/", figpath::String = "outputfig/",  levels::Array = [0, 0, 0], cmap::String = "jet", cbar::Bool = true, plotdim::String = "2D")
    for frame in frames
        outputfig(frame, varname = varname, filepath = filepath, figpath = figpath, levels = levels, cmap = cmap, cbar = cbar, plotdim = plotdim)
    end
end

"""
从内存中读取远比硬盘快。但也需要写一个从硬盘读取并绘图的函数。
"""
function outputfig!(f::Fluid; varname::String = "rho", axis::Int = 1, frame::Int = 1, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,1,1])
    if plotdim in ("2D", "2D-x", "2D-y")
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        y = [f.d[2]*(i-0.5-f.ng)+f.point1[2]  for  i=1:f.nmesh[2]+f.ng*2]
        if varname in ("rho","p","e")
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k], Symbol(varname))
            end
        elseif varname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end  
        else
            error("undef varname")          
        end  
        ioff()
        if plotdim == "2D-x"
            output_plot_func(x, data[:, floor(Int, size(data, 2)/2)], levels=levels)
        elseif plotdim == "2D-y"
            output_plot_func(y, data[floor(Int, size(data, 1)/2), :], levels=levels)  
        else        
            output_plot_func(x, y, data, levels =levels)
        end
    elseif plotdim == "1D"
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        data = Array{Float64}(undef, length(x))
        for k in eachindex(f.cells)
            data[k] = getfield(f.cells[k], Symbol(varname))
        end
        ioff()
        output_plot_func(x, data, levels=levels)        
    else
        error("undef plotdim")
    end
    savefig(figpath*varname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end

function pid_to_index(pid::Int; nzone::Array{Int} = [1], dist::Array{Int} = [1])
    nsize = nzone ./ dist
    dim = length(nzone)
    I = zeros(Int, dim)
    if dim >= 1
        I[1] = pid%dist[1]
    end
    if dim >= 2
        I[2] = Int( cld(pid%prod(dist[1:2]), dist[1]) )
    end
    if dim >= 3
        I[3] = Int( cld(pid, prod(dist[1:2])) )
    end
    return Tuple([Int(nsize[k]*(I[k]-1)+1):Int(nsize[k]*I[k]) for k in 1:dim])    
end

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

function read_frame(n::Int; varname::String = "rho", filepath::String = "../outputdata/")
    var = readdlm(filepath*"fluid_"*varname*"_"*string(n+FRAME_BASE)*".txt")
    return var
end

function read_mesh(dim::Int, n::Int; varname::String = "rho", filepath::String = "../outputdata/")
    return Tuple([readdlm(filepath*"fluid_x"*string(k)*".txt") for k in 1:dim])
end

function reflective_cell!(c::GCell, axis::Int)
    rc = copy!(c)
    rc.w[1+axis] *= -1
    rc.wb[1+axis] *= -1
    rc.u[axis] *= -1
    return rc
end

function reynolds_number(rho::Float64, constants::Dict)
    return rho * constants["U0"] * constants["L0"] / constants["mu"]
end

function set_boundaries!(f::Fluid, boundaries::Array{String})
    f.boundaries = boundaries
    update_boundaries!(f)
end

function showdist!(a::DArray)
    println("PID   |  localindices")
    for pid in workers()
        ind = @fetchfrom pid localindices(a)
        for i in ind
            if length(i) == 0
                ind = "-"
                break
            end
        end
        println(pid,"        ",ind)
    end
end

function showfield!(c::AbstractArray, str::String, inds::UnitRange{Int64}...; axis::Int = 0)
    s = Symbol(str)
    # 如果没有输入inds则输出全部。
    if length(inds) == 0
        inds = Tuple([1:size(c,i) for i = 1:length(size(c))])
    end
    a = Array{Any}(undef, size(view(c, inds...))) 
    firstinds = CartesianIndex(Tuple([inds[i][1] - 1 for i = 1:length(inds)]))
    if axis == 0
        for k in CartesianIndices(a)
            a[k] = getfield(c[k + firstinds], s)
        end
    else
        for k in CartesianIndices(a)
            ak = getfield(c[k + firstinds], s)
            a[k] = ak[axis]
        end
    end
    display(a)
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

function solve!(f::Fluid; CFL::Float64 = 0.3, maxtime::Float64 = 1, maxframe::Int = 1, cutframe::Int = 1, varname::String = "rho", axis::Int = 1, filepath::String = "../outputdata/", draw::Bool = false, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,0,0])
    time = 0.
    frame = 0
    dt = 0.
    println("|||||||||||||||||||||||||||||||||||||||||||")
    @printf "Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100

    if OUTPUTDATA
        output!(Int(frame/cutframe), time, filepath = filepath)
        output!(f, varname = varname, axis = axis, frame = Int(frame/cutframe), filepath = filepath)
        if draw
            outputfig!(f, frame = Int(frame/cutframe), varname = varname, axis = axis, figpath = figpath, plotdim = plotdim, levels=levels)
        end
    end
    while time < maxtime && frame < maxframe
        dt = time_step!(f, CFL = CFL)
        fvm_advance!(f, dt)
        time += dt
        frame += 1
        print("Current frame = ",frame,"\r")
        if frame%cutframe == 0
            println("|||||||||||||||||||||||||||||||||||||||||||")
            @printf "Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100
            if OUTPUTDATA
                output!(Int(frame/cutframe), time, filepath = filepath)
                output!(f, varname = varname, axis = axis, frame = Int(frame/cutframe), filepath = filepath)
                if draw
                    outputfig!(f, frame = Int(frame/cutframe), varname = varname, axis = axis, figpath = figpath, plotdim = plotdim,levels=levels)
                end                
            end
        end
    end
end

function states2w(;rho::Float64=1.0,u::Array{Float64}=[0.], e::Float64=1.0)
    return [rho; rho * u; rho * (e+0.5*MK.norm2(u))]
end

function temperature(c::GCell, gamma::Float64)
    s = sound_speed(rho = c.rho, p = c.p, gamma = gamma)
    if isnan(s)
        println("-- temperature --")
        display(c)
        error("T error")
    end
    return s^2
end

"""
像time_step!()这样仅读取的函数，则无需远程，可直接读取DArray。
"""
function time_step!(f::Fluid; CFL::Float64 = 0.1)
    smax = 0.
    @sync for pid in workers()
        smax = max(smax, @fetchfrom pid begin
            smaxtmp = 0.
            for c in localpart(f.cells)
                if c.rho == 0.
                    s = 0.
                else
                    s = sound_speed(rho = c.rho, p = c.p, gamma = f.constants["gamma"])
                end
                smaxtmp = maximum([smaxtmp; s .+ map(abs, c.u)])
            end
            smaxtmp
        end)
    end
    dt = minimum(f.d) / smax * CFL
    if dt < TOL_STEP
        error( "Too small time step!")
    end
    return dt
end

"""
一个合理的思路是：
1. 把要修改的子数组的索引modrange发送到pid分区。

2. 在pid分区上操作。首先确定要修改的子数组索引modrange与pid分区索引pidrange的交集intrange。

3. 对于不需要修改的子数组，直接从本地拉取到pid分区。

4. 其余不需要修改的参数也从本地拉取即可。

5. 在pid分区调用函数计算，最后把结果拉取到本地。
"""
function update_boundaries!(f::Fluid)
    nsize = (f.nmesh .+ f.ng*2) ./f.dist
    for n in nsize
        if n < f.ng * 2
            error("Too small nsize")
        end
    end
    @sync for pid in workers() 
        @spawnat pid begin
            pid_range = localindices(f.cells)
            bias = [r[1] for r in pid_range] .- 1
            for axis in 1:f.dim
                for side in 1:2
                    for i in (f.ng+f.nmesh[axis])*(side-1) .+ (1:f.ng)
                        int_range = map(intersect, pid_range, Tuple([k == axis ? (i:i) : (1:size(f.cells, k)) for k in 1:f.dim]))
                        if length_of_tuples(int_range) > 0  # 如果intrange索引的“覆盖面积”是0,那么就不做计算。
                            ref_i = 2*(f.ng+f.nmesh[axis]*(side-1))+1-i
                            ref_range = Tuple([k == axis ? (ref_i:ref_i) : int_range[k] for k = 1:f.dim])
                            local_int_range = Tuple([int_range[k] .- bias[k] for k in 1:f.dim])
                            ghost_boundary_cell!(f.cells[ref_range...], localpart(f.cells)[local_int_range...], f.boundaries[axis, side], axis = axis)   
                            # println(("side",side,"pid_range = ",pid_range,"bias = ",bias,"int_range = ",int_range,"local_int_range",local_int_range,"ref_range = ",ref_range))                                                
                        end
                    end
                end
            end
        end
    end
    # showfield!(f.cells, "rho")
end

function update_cells!(f::Fluid, rk::Int, dt::Float64)
    coeff = RK_COEFF[:, rk]
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                if MK.betweeneq(c.x, f.point1, f.point2)
            
                    @. c.w = coeff[1] * c.w + coeff[2] * c.wb + coeff[3] * c.rhs * dt


                    # correct_cell_w!(c, f.constants["gamma"])
                    w2states!(c, f.constants["gamma"])  

                end
            end
        end
    end
end

function update_rhs!(f::Fluid)
    if f.dim == 1
        @sync for pid in workers()
            @spawnat pid begin
                ws = Array{Float64,2}(undef, 2+f.dim, 5)
    
                for c in localpart(f.cells) # 千万别漏写localpart 
                    if MK.betweeneq(c.x, f.point1, f.point2)
                        c.rhs = zeros(Float64, 2+f.dim)
                        ci = copy(c.i)
    
                        for j in 1:5
                            sw = f.cells[ci[1]+j-3].w
                            for i = 1:2+f.dim
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 1) / f.d[1] 
                                         
                    end
                end
            end
        end
    elseif f.dim == 2
        @sync for pid in workers()
            @spawnat pid begin
                ws = Array{Float64,2}(undef, 2+f.dim, 5)
    
                for c in localpart(f.cells) # 千万别漏写localpart 
                    if MK.betweeneq(c.x, f.point1, f.point2)
                        c.rhs = zeros(Float64, 2+f.dim)
                        ci = copy(c.i)
    
                        for j in 1:5
                            sw = f.cells[ci[1]+j-3,ci[2]].w
                            for i = 1:4
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 1) / f.d[1] 
                        
    
                        for j in 1:5
                            sw = f.cells[ci[1],ci[2]+j-3].w
                            for i = 1:4
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 2) / f.d[2]                    
                    end
                end
            end
        end
    else
        error("undef dim")
    end
end 

# """
# 这样写会慢30%
# """
# function testupdate_rhs!(f::Fluid)
#     @sync for pid in workers()
#         @spawnat pid begin
#             ws = Array{Float64,2}(undef, 2+f.dim, 5)
#             for c in localpart(f.cells) # 千万别漏写localpart 
#                 if MK.betweeneq(c.x, f.point1, f.point2)     
#                     c.rhs = zeros(Float64, 2+f.dim)
#                     for axis = 1:f.dim
#                         i = copy(c.i)
                         
#                         for k = 1:5
#                             i[axis] = c.i[axis] + k - 3
#                             ws[:,k] = f.cells[Tuple(i)...].w
#                         end 
                        
#                         c.rhs += (get_flux!(f.reconst_scheme, ws[:,1:4], f.constants, axis = axis) - get_flux!(f.reconst_scheme, ws[:,2:5], f.constants, axis = axis)) / f.d[axis]
#                     end                       
                    
#                 end
#             end
#         end
#     end
# end 

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

function w2states!(c::GCell, gamma::Float64)
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

###########
end