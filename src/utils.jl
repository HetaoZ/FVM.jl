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

function get_point_LD!(x::Vector{Float64}, f::Fluid)
    return [ceil(Int, (x[k]-f.point1[k])/f.d[k]-0.5+f.ng) for k = 1:f.dim]
end

function get_point_RU!(x::Vector{Float64},f::Fluid)
    iLD = get_point_LD!(x, f)
    return iLD .+ 1
end

function get_point_located_cell!(x::Vector{Float64}, f::Fluid)
    i = get_point_LD!(x, f)
    ip = copy(i)
    for k in 1:f.dim
        if x[k] > (i[k] - f.ng) * f.d[k] + f.point1[k]
            ip[k] += 1
        end
    end  
    return ip
end