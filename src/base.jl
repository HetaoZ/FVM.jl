

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
    ndiv::Int # 分区总数 = prod(dist)
    # 单元
    cells::DArray{Cell,N,A} where N where A
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
    reconst_scheme::String
    flux_scheme::String
end

"""
创建一个流场。
"""
function Fluid(dim::Int; point1::Array = [0.], point2::Array = [1.], nmesh::Array = [1], ng::Int = 1, dist::Array = [1], constants::Dict = Dict())
    # if (nmesh .+ ng*2) .% dist != zeros(Int, length(nmesh))
    #     error("(nmesh .+ ng * 2) must be divisible by dist.")
    # end
    # println(typeof(Cell(1)))
    Fluid(dim, point1, point2, nmesh, (point2 - point1) ./ nmesh, ng, dist, prod(dist),
    # 分布式矩阵
    distribute([Cell(dim)]), 
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
    false, 0., false, false,
    "MUSCL", "AUSM"
    )
end



