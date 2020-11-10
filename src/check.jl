function check_cell_status!(c::Cell)
    nan = false
    if isnan(sum(c.w)) || isnan(sum(c.rhs)) || isnan(sum(c.wb))
        nan = true
    end

    if nan
        println("-- c is NaN --")
        display(c)
        error("NaN error")
    end
end

function check_all_cells_status!(cells::DArray)
    for c in cells
        check_cell_status!(c)
    end
end

function check_all_cells_status!(f::Fluid)
    for c in f.cells
        check_cell_status!(c)
    end
end

function check_all_cells_particles_status!(f::Fluid)
    check_all_cells_status!(f.cells)
    check_all_particles_status!(f.particles)
end

function check_all_particles_status!(particles::Array{Particle})
    for p in particles
        check_particle_status!(p)
    end
end

function check_all_particles_status!(f::Fluid)
    for p in f.particles
        check_particle_status!(p)
    end
end

function check_conservativity!(f::Fluid)
    @warn "This may take much time!"
    mass = 0.
    for c in f.cells
        if MK.betweeneq(c.x, f.point1, f.point2)
            mass += c.rho * MK.product(f.d) 
        end
    end
    if !f.total_is_summed
        f.total_mass = mass
        f.total_is_summed = true
    end
    println("-- [ mass = ", mass, " | δ = ", (mass-f.total_mass)/f.total_mass, " ] --")    
end

function check_mass!(f::Fluid; point1::Array{Float64} = f.point1, point2::Array{Float64} = f.point2)
    mass = 0.
    for c in f.cells
        if MK.betweeneq(c.x, point1, point2)
            mass += c.rho * MK.product(f.d) 
        end
    end
    if !f.total_is_summed
        f.total_mass = mass
        f.total_is_summed = true
    end
    
    return mass
end

function check_mass!(particles::Array{Particle})
    mass = 0.
    for p in particles
        mass += p.m
    end
    return mass
end

function check_particle_status!(p::Particle)
    nan = false
    if isnan(sum(p.m)) || isnan(sum(p.mu)) || isnan(sum(p.mE))
        nan = true
    end
    if isnan(sum(p.x)) || isnan(sum(p.dx)) || isnan(sum(p.force_to_boundary))
        nan = true
    end
    if isnan(sum(p.boundary_id)) || isnan(sum(p.boundary_ratio))
        nan = true
    end

    if nan
        println("-- p is NaN --")
        display(p)
        error("NaN error")
    end    
end

"""
在密度相对较低的单元，数值伪振荡可能造成负质量或负压力。对其进行保正性修正，会造成轻微不守恒，总质量会轻微增加。
"""
function correct_cell_w!(c::Cell, gamma)
    if sign.([c.w[1], c.w[end], pressure(c.w, gamma)]) == [-1.0, -1.0, -1.0]
        # 允许负质量单元存在。
    elseif c.w[1] < 0 || c.w[end] < 0 || pressure(c.w, gamma) < 0
        # println("-- correct_cell_w: 1 --")
        # println("c.i = ", c.i)

        c.w = zeros(Float64, length(c.w))
    end
end
