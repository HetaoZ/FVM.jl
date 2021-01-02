function get_point_LD!(x::Vector{Float64}, f::Fluid)
    id = ones(Int, 3)
    id[1:f.realdim] = [ceil(Int, (x[k]-f.point1[k])/f.d[k]-0.5+f.ng) for k = 1:f.realdim]
    return id
end

function get_point_RU!(x::Vector{Float64},f::Fluid)
    iLD = get_point_LD!(x, f)
    return iLD .+ 1
end

function get_point_located_cell!(x::Vector{Float64}, f::Fluid)
    i = get_point_LD!(x, f)
    ip = copy(i)
    for k in 1:f.realdim
        if x[k] > (i[k] - f.ng) * f.d[k] + f.point1[k]
            ip[k] += 1
        end
    end  
    return ip
end

# function check_mass!(f::Fluid)
#     id = [f.ng+1:f.ng+f.nmesh[1], f.ng+1:f.ng+f.nmesh[2], f.ng+1:f.ng+f.nmesh[3]]
#     if f.realdim < 3
#         id[3] = 1:1
#     end
#     if f.realdim < 2
#         id[2] = 1:1
#     end

#     return sum(f.rho[id[1], id[2], id[3]]) * prod(f.d[1:f.realdim])
# end

function check_marked_mass!(f::Fluid, marker::Int)
    
    V = prod(f.d[1:f.realdim])

    mass = @sync @distributed (+) for id in CartesianIndices(f.rho)
        local_mass = 0.
        if f.marker[id] == marker
            if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
                local_mass += f.rho[id] * V
            end
        end
        local_mass
    end

    return mass
end

function check_marked_mass!(f::Fluid, markers::Vector{Int})
    
    V = prod(f.d[1:f.realdim])

    mass = @sync @distributed (+) for id in CartesianIndices(f.rho)
        local_mass = 0.
        if f.marker[id] in markers
            if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
                local_mass += f.rho[id] * V
            end
        end
        local_mass
    end

    return mass
end

function check_mass!(f::Fluid)
    
    V = prod(f.d[1:f.realdim])

    mass = @sync @distributed (+) for id in CartesianIndices(f.rho)
        local_mass = 0.
        if f.marker[id] != 0
            if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
                local_mass += f.rho[id] * V
            end
        end
        local_mass
    end

    return mass
end

function correct_cell_w(w, gamma, rho0, u0, e0)

    rho, u, e, p = w_to_status(w, gamma)
    if rho < -1e-14
        println(w)
        println((rho,u,e,p))
        error("rho < 0 ")
        return zeros(Float64, 5)
    end
    
    if isnan(rho+e+p+sum(u))
        println(w)
        println((rho,u,e,p))
        error("isnan(rho+e+p+sum(u)) ")
        return status_to_w(rho0, u0, e0)
    end
    if e < 0.
        # println(w)
        # println((rho,u,e,p))
        # error("e < 0. ")
        return status_to_w(rho, u, e0*1e-7)
    end
    if e > 1e7*e0
        println(w)
        println((rho,u,e,p))
        error("e > 1e7*e0 ")
        return status_to_w(rho0, u0, e0)
    end
    if pressure(w, gamma) < 0.
        # println(w)
        # println((rho,u,e,p))
        # error("pressure(w, gamma) < 0. ")
        return status_to_w(rho, u, e0*1e-7)
    end
    return w
end