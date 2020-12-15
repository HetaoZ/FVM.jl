function length_of_tuples(r::Tuple)
    s = 1
    for k in r
        s *= length(k)
    end
    return s
end

function index_to_pid(k::Union{Array{Int}, CartesianIndex}; nzone::Array{Int} = [1], dist::Array{Int} = [1])
    K = cld.(Tuple(k), nzone ./ dist)
    pid = workers()[Int(K[1] + [K[i] - 1 for i in 2:length(K)]' * [prod(dist[1:i]) for i in 1:length(K)-1])]
    return pid
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

