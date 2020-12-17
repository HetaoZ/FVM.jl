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

function check_mass!(f::Fluid)
    id = [f.ng+1:f.ng+f.nmesh[1], f.ng+1:f.ng+f.nmesh[2], f.ng+1:f.ng+f.nmesh[3]]
    if f.realdim < 3
        id[3] = 1:1
    end
    if f.realdim < 2
        id[2] = 1:1
    end

    return sum(f.rho[id[1], id[2], id[3]]) * prod(f.d[1:f.realdim])
end