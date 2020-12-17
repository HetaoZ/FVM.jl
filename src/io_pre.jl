"""
向sdir方向运动的激波，波前是区域1，波后是区域2，激波本身记为s，激波马赫数Ms=激波的地面速度us的绝对值除以波前声速c1。

算例：Monasse(2012) Example 2

mshock(1.01e5,1.2,0,1.21,1.4,-1) = (1.630779226806516, 155686.45, -109.71829086299334)
"""
function after_shock(p1,ρ1,u1,Ms,γ,sdir)
    @assert Ms >= 1
    c1 = sound_speed(ρ1, p1, γ)
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

function fill_fluid!(f::Fluid, rho, u, p)
    f.para["rho0"] = rho
    f.para["u0"] = u
    f.para["p0"] = p
    e = pressure_to_e(rho, p, f.para["gamma"])
    f.para["e0"] = e
    fill_fluid!(f, f.point1, f.point2, rho, u, p)
    f.para["background is filled"] = true  
end

function fill_fluid!(f::Fluid, point1::Array, point2::Array, rho, u, p)
    e = pressure_to_e(rho, p, f.para["gamma"])
    w = status_to_w(rho, u, e)
    @sync @distributed for id in CartesianIndices(f.rho)
        if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], point1, point2)
            f.rho[id] = rho
            f.u[:,id] = u
            f.e[id] = e
            f.p[id] = p
            f.w[:,id] = w
        end
    end 
end

function clear_fluid_in_box!(f, point1, point2)
    fill_fluid!(f, point1, point2, 0., [0., 0., 0.], 0.)
end

function set_bounds!(f, bound_types...)
    f.boundx = bound_types[1]
    if f.realdim > 1
        f.boundy = bound_types[2]
    end
    if f.realdim > 2
        f.boundz = bound_types[3]
    end
    update_bounds!(f)
end

"""
Use 'command | tee -a log' in Linux Terminal to make a log.
"""
function review(f::Fluid)
    println("-- short summary of Fluid --")
    println("# parameters")
    for k in keys(f.para)
        println("  ",k," : ",f.para[k])
    end
    println("# mesh")
    println("  real dimension : ", f.realdim)
    println("  domain : ", Tuple(f.point1)," -> ", Tuple(f.point2))
    println("  number of cells : ", length(f.rho))
    println("  size of cells : ", size(f.rho))
    println("# distribution")
    println("  number of workers : ", nworkers())
    println("# boundaries")
    println("  axis 1 : ", f.boundx)
    println("  axis 2 : ", f.boundy)
    println("  axis 3 : ", f.boundz)
    
    println("# physical status")
    println("  rho ∈  ", [minimum(f.rho), maximum(f.rho)])
    println("  p ∈  ", [minimum(f.p), maximum(f.p)])
end