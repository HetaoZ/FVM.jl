
# using Distributed
# addprocs(1)

println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere using FVM

@everywhere using DistributedArrays

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(realdim = 2, 
            point1 = [-10e-3, 0, 0], 
            point2 = [50e-3, 65e-3, 0], 
            nmesh = Int[60, 65, 1] .* 1, 
            ng = 2, 
            dist = [nworkers(), 1, 1])  
f.para["gamma"] = 1.4
f.para["viscosity"] = false
f.para["flux scheme"] = "LF" # "AUSM" or "LF"

rho0, u0, p0 = 1.0, [0., 0., 0.], 1.0
FVM.fill_fluid!(f, rho0, u0, p0)

rho2, p2, u2 = FVM.after_shock(p0, rho0, u0[1], 1.21, f.para["gamma"], 1)

FVM.fill_fluid!(f, [-10e-3, 0., 0.], [0, 65e-3, 0.], rho2, [u2, 0., 0.], p2)

FVM.clear_fluid_in_box!(f, [0, 0, 0], [5e-3, 50e-3, 0])

FVM.set_bounds!(f, ["free", "refl"], ["refl", "refl"])

# review(f)

# # --------------------------------
# # solve
# # --------------------------------
frame = 0
time = 0
N = 1000000

while frame < 1 && time < 1e-3
    global m, frame, time
    dt = FVM.time_step!(f, CFL=0.05)
    FVM.advance!(f, dt)
    frame += 1
    time += dt
    if frame%1 == 0
        FVM.save_to_vtk(f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
        println(frame,"   ",dt, "   ",time)
    end
end