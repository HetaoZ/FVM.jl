mutable struct Particle
    dim::Int
    m::Float64 # 质量
    u::Array{Float64} # 动量
    ek::Float64 # 动能
    e::Float64 # 内能
    E::Float64 # 总能
    x::Array{Float64} # 坐标
    dx::Array{Float64} # 运动矢量
    V::Float64 # volume
    r::Float64 # radius
    force_to_boundary::Array{Float64} # 对边界的力
    boundary_id::Int # 力所在的边界编号
    boundary_ratio::Float64 # 力的作用点在boundary上的比例（左/全）
end
Particle(dim::Int) = Particle(dim, 0., zeros(Float64,dim), 0., 0., 0., zeros(Float64,dim), zeros(Float64,dim), 0., 0., zeros(Float64,dim), 0, 0.)

mutable struct Cell 
    # 物理量
    rho::Float64
    u::Array{Float64,1}
    e::Float64
    p::Float64
    # 迭代量
    w::Array{Float64,1}
    wb::Array{Float64,1}
    rhs::Array{Float64,1}
    # 几何量
    i::Array{Int,1} # 网格编号
    x::Array{Float64,1}
end

function Cell(dim::Int; rho::T = 1.0, u::Array{T} = zeros(Float64,dim), p::T = 1., gamma::Float64 = 1.4) where T <: Real
    Cell(
    # 物理量
    rho, u, pressure_to_e(rho=rho, p=p, gamma=gamma), p,
    # 迭代量
    states2w(rho=rho,u=u,e=pressure_to_e(rho=rho, p=p, gamma=gamma)),
    states2w(rho=rho,u=u,e=pressure_to_e(rho=rho, p=p, gamma=gamma)), 
    zeros(Float64,2+dim), 
    # 几何量
    zeros(Int,dim), zeros(Float64,dim))  
end

"""
流场
"""
mutable struct Fluid
    # 维数
    dim::Int
    # 几何量
    point1::Array{T} where T <: Real 
    point2::Array{T} where T <: Real 
    nmesh::Array{Int} # 各维度网格数
    d::Array{T} where T <: Real# 各维度步长
    ng::Int 
    dist::Array{Int} # 各维度分区数向量
    ndiv::Int # 分区总数 = MK.product(dist)
    # 单元
    cells::DArray{Cell,N,A} where N where A
    # 虚拟网格单元
    ghost_cells::Dict{Tuple, Cell}
    # 粒子
    particles::Array{Particle}
    # 边界类型
    boundaries::Array{String}
    # 可以像C++一样放一个成员函数在这里，但没必要。
    # solver::Function
    # 常数
    constants::Dict
    # 辅助变量
    background_is_filled::Bool
    total_mass::Float64
    total_is_summed::Bool
    consider_vis_item::Bool
    exclude_particles::Bool
    reconst_scheme::String
    flux_scheme::String
end

"""
创建一个流场。
"""
function Fluid(dim::Int; point1::Array = [0.], point2::Array = [1.], nmesh::Array = [1], ng::Int = 1, dist::Array = [1], constants::Dict = Dict())
    if (nmesh .+ ng*2) .% dist != zeros(Int, length(nmesh))
        error("(nmesh .+ ng * 2) must be divisible by dist.")
    end
    # println(typeof(Cell(1)))
    Fluid(dim, point1, point2, nmesh, (point2 - point1) ./ nmesh, ng, dist, MK.product(dist),
    # 分布式矩阵
    distribute([Cell(dim)]), 
    # 虚拟网格单元
    Dict(),
    # 粒子
    Particle[],
    # 边界
    Array{String}(undef, dim, 2),
    # 常数
    Dict("gamma"=>1.4, # 多方指数
    "mu"=>1.e-6, # 动力粘性系数
    "Pr"=>1.0, # 热传导系数
    "L0"=>1, # 特征长度
    "U0"=>1 # 特征速度
    ),
    # 辅助变量
    false, 0., false, false, true,
    "MUSCL", "AUSM"
    )
end

function length_of_tuples(r::Tuple)
    s = 1
    for k in r
        s *= length(k)
    end
    return s
end

function index_to_pid(k::Union{Array{Int}, CartesianIndex}; nzone::Array{Int} = [1], dist::Array{Int} = [1])
    K = cld.(Tuple(k), nzone ./ dist)
    pid = workers()[Int(K[1] + [K[i] - 1 for i in 2:length(K)]' * [MK.product(dist[1:i]) for i in 1:length(K)-1])]
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
        I[2] = Int( cld(pid%MK.product(dist[1:2]), dist[1]) )
    end
    if dim >= 3
        I[3] = Int( cld(pid, MK.product(dist[1:2])) )
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