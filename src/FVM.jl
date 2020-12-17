
module FVM
using DelimitedFiles, Printf, Distributed, LinearAlgebra, Statistics, WriteVTK
using SharedArrays
using DUtils
using MathKits; const MK = MathKits
export Fluid

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
const BOUND_TYPE = Dict("free"=>1.0, "refl"=>-1.0)

include("base.jl")
include("utils.jl")
include("physical.jl")
include("io_pre.jl")
include("boundary.jl")
include("flux.jl")
include("solver.jl")
include("copy.jl")
include("io_post_vtk.jl")

###
end