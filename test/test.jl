

nw = 4

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())
try
    @everywhere using .FVM
catch
    @everywhere include("src/FVM.jl")
    @everywhere using .FVM
end
@everywhere using DistributedArrays, MathKits

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-10e-3,-10e-3], 
            point2 = [20e-3, 65e-3], 
            nmesh = Int[30, 75] .* 4, 
            ng = 4, 
            dist = [nw, 1]
            )  
f.constants["gamma"] = 1.4
f.consider_vis_item = false

c1 = Cell(2, rho = 1., u = [0., 0.], p = 1.)
fill_fluid!(f, c1)

rho2, p2, u2 = after_shock(c1.p,c1.rho,c1.u[1],1.21,f.constants["gamma"],1)

c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
fill_fluid!(f, c2, [-10e-3, 0.], [0., 65e-3])

@everywhere function clear_cell!(c::Cell)
    r = 0.0
    c.rho = c.rho*r
    c.u = zeros(Float64, size(c.u))
    c.e = c.e
    c.p = c.p*r
    c.w = [c.rho, 0, 0, c.rho*c.e]
end
function clear_fluid_in_box!(f, point1, point2)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                if MathKits.between(c.x, point1, point2)
                    clear_cell!(c) 
                end
            end
        end
    end
end
clear_fluid_in_box!(f, [0, 0], [5e-3, 50e-3])

FVM.set_bounds!(f, ["free" "refl"; "refl" "refl"])

f.exclude_particles = false
f.flux_scheme = "LF"


# --------------------------------
# solve
# --------------------------------
frame = 0
time = 0
N = 1000000

while frame < 100 && time < 10.e-3
    global m, frame, time
    dt = FVM.time_step!(f, CFL=0.05)
    FVM.advance!(f, dt)
    frame += 1
    time += dt
    save_to_vtk(f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
    println(frame,"   ",dt, "   ",time)
end