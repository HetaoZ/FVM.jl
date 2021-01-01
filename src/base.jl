"""
3D fluid definition
"""
mutable struct Fluid
    realdim::Int

    point1::Vector{Real}
    point2::Vector{Real}
    nmesh::Vector{Int}
    d::Vector{Real}
    ng::Int

    NX::Int
    NY::Int
    NZ::Int

    x::Vector{Real}
    y::Vector{Real}
    z::Vector{Real}

    rho::SharedArray{Float64,3}
    u::SharedArray{Float64,4}
    e::SharedArray{Float64,3}
    p::SharedArray{Float64,3}

    w::SharedArray{Float64,4}
    wb::SharedArray{Float64,4}
    rhs::SharedArray{Float64,4}

    marker::SharedArray{Int8,3}

    boundx::Vector{String}
    boundy::Vector{String}
    boundz::Vector{String}

    para::Dict
end

function Fluid(; realdim::Int = 3, point1::Array = [0, 0, 0], point2::Array = [1, 1, 1], nmesh::Vector{Int} = [1, 1, 1], ng::Int = 2, para::Dict = Dict())

    default_para = Dict("rho0"=>0.,
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
    "flux scheme"=>"AUSM",
    "fill forever"=>[]
    )

    for key in keys(para)
        default_para[key] = para[key]
    end

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

    return Fluid(realdim, point1, point2, nmesh, d, ng,

    NX, NY, NZ,
    x, y, z,

    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (3, NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 

    SharedArray(zeros(Float64, (5, NX, NY, NZ))),
    SharedArray(zeros(Float64, (5, NX, NY, NZ))),
    SharedArray(zeros(Float64, (5, NX, NY, NZ))),

    SharedArray(ones(Int8, (NX, NY, NZ))),

    ["none", "none"],
    ["none", "none"],
    ["none", "none"],

    default_para)
end