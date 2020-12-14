"""
3D fluid definition
"""
mutable struct Fluid
    realdim::Int

    point1::Array 
    point2::Array
    nmesh::Vector{Int}
    d::Vector
    ng::Int

    dist::Array{Int} 
    ndiv::Int

    NX::Int
    NY::Int
    NZ::Int

    x::Vector
    y::Vector
    z::Vector 

    rho::DArray{Float64}
    u::DArray
    e::DArray{Float64}
    p::DArray{Float64}

    w::DArray
    wb::DArray
    rhs::DArray

    mark::DArray{Int8}
    target_id::DArray

    boundx::Array{String}
    boundy::Array{String}
    boundz::Array{String}

    para::Dict
end

function Fluid(; realdim::Int = 3, point1::Array = [0, 0, 0], point2::Array = [1, 1, 1], nmesh::Vector{Int} = [1, 1, 1], ng::Int = 2, dist::Array = [1, 1, 1], para::Dict = Dict("rho0"=>0.,
    "u0"=>[0., 0., 0.],
    "e0"=>0.,
    "p0"=>0.,
    "gamma"=>1.4, # 多方指数
    "mu"=>1.e-6, # 动力粘性系数
    "Pr"=>1.0, # 热传导系数
    "L0"=>1, # 特征长度
    "U0"=>1, # 特征速度
    "background is filled"=>false, 
    "total mass"=>0., 
    "total is summed"=>false, 
    "viscosity"=>false,
    "reconst scheme"=>"MUSCL", 
    "flux scheme"=>"AUSM"
    ))

    d = (point2 - point1) ./ nmesh

    NX = nmesh[1]+2*ng
    x = [point1[1] + d[1] * (i-0.5-ng) for i = 1:NX]
    
    if realdim > 1
        NY = nmesh[2]+2*ng
        y = [point1[2] + d[2] * (i-0.5-ng) for i = 1:NY]
    else
        NY = 1
        y = [point1[2]]
    end
    if realdim > 2
        NZ = nmesh[3]+2*ng
        z = [point1[3] + d[3] * (i-0.5-ng) for i = 1:NZ]
    else
        NZ = 1
        z = [point1[3]]
    end

    return Fluid(realdim, point1, point2, nmesh, d, ng, dist, prod(dist),

    NX, NY, NZ,
    x, y, z,

    distribute(fill(para["rho0"], (NX, NY, NZ)), procs = workers(), dist = dist), 
    distribute(fill(para["u0"], (NX, NY, NZ)), procs = workers(), dist = dist), 
    distribute(fill(para["e0"], (NX, NY, NZ)), procs = workers(), dist = dist), 
    distribute(fill(para["p0"], (NX, NY, NZ)), procs = workers(), dist = dist), 

    distribute(fill(zeros(Float64, 5), (NX, NY, NZ)), procs = workers(), dist = dist),
    distribute(fill(zeros(Float64, 5), (NX, NY, NZ)), procs = workers(), dist = dist),
    distribute(fill(zeros(Float64, 5), (NX, NY, NZ)), procs = workers(), dist = dist),

    distribute(fill(Int8(1), (NX, NY, NZ)), procs = workers(), dist = dist),
    distribute(fill(zeros(Int, 3), (NX, NY, NZ)), procs = workers(), dist = dist),

    ["none", "none"],
    ["none", "none"],
    ["none", "none"],

    para)
end