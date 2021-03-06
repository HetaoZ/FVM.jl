
# using Distributed
# addprocs(1)

println("Opened ", nworkers()," process(es) of PID ", workers())

try 
    @everywhere using .FVM
catch
    @everywhere include("src/FVM.jl")
    @everywhere using .FVM
end

@everywhere using DistributedArrays

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
# f = Fluid(realdim = 2, 
#             point1 = [-10e-3, 0, 0], 
#             point2 = [50e-3, 65e-3, 0], 
#             nmesh = Int[60*4, 65*4, 1], 
#             ng = 2)
# f.para["gamma"] = 1.4
# f.para["viscosity"] = false
# f.para["flux scheme"] = "LF" # "AUSM" or "LF"

# rho0, u0, p0 = 1.0, [0., 0., 0.], 1.0
# FVM.fill_fluid!(f, rho0, u0, p0)

# rho2, p2, u2 = FVM.after_shock(p0, rho0, u0[1], 1.21, f.para["gamma"], 1)

# FVM.fill_fluid!(f, [-10e-3, 0., 0.], [0, 65e-3, 0.], rho2, [u2, 0., 0.], p2)

# # FVM.clear_fluid_in_box!(f, [0, 0, 0], [1e-3, 1e-3, 0])

# FVM.set_bounds!(f, ["free", "free"], ["refl", "refl"])

# FVM.review(f)

f = Fluid(
    realdim = 2, 
    point1 = [-1e-3, 0, 0], 
    point2 = [3e-3, 2e-3, 0], 
    nmesh = Int[4*10, 2*10, 1], 
    ng = 2
)  
f.para["gamma"] = 1.4
f.para["viscosity"] = false
f.para["flux scheme"] = "LF" # "AUSM" or "LF"

rho0, u0, p0 = 1.0, [0., 0., 0.], 1.0
FVM.fill_fluid!(f, rho0, u0, p0)

# rho2, p2, u2 = FVM.after_shock(p0, rho0, u0[1], 1.21, f.para["gamma"], 1)

FVM.fill_fluid!(f, [0., 0., 0.], [1e-3, 1e-3, 0.], 0., [0., 0., 0.], 0.)

FVM.set_bounds!(f, ["free", "refl"], ["refl", "refl"])

# # --------------------------------
# # solve
# # --------------------------------
function test(f, maxframe)
    frame = 0
    time = 0
    N = 1000000

    println("mass = ", FVM.check_mass!(f))
    FVM.save_to_vtk(f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
    println(frame,"   ",0, "   ",time)

    while frame < maxframe 
        
        dt = FVM.time_step!(f, CFL=0.05)
        FVM.advance!(f, dt)
        println("mass = ", FVM.check_mass!(f))
        frame += 1
        time += dt
        if frame%1 == 0
            FVM.save_to_vtk(f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
            println(frame,"   ",dt, "   ",time)
        end
    end

end

test(f, 5)