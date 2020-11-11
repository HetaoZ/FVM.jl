
module FVM
using PyPlot 
using DelimitedFiles, Printf, Distributed, DistributedArrays, LinearAlgebra, Statistics
using MathKits
const MK = MathKits
export Fluid, Cell, fill_fluid!, fvm_set_bounds!, after_shock, fvm_solve!

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


include("base.jl")
include("io_pre.jl")
include("state.jl")
include("boundary.jl")
include("flux.jl")
include("source.jl")
include("solver.jl")
include("deepcopy.jl")
include("check.jl")
include("io_post.jl")


###########
end