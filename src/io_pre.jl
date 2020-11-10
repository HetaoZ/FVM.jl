"""
向sdir方向运动的激波，波前是区域1，波后是区域2，激波本身记为s，激波马赫数Ms=激波的地面速度us的绝对值除以波前声速c1。

算例：Monasse(2012) Example 2

mshock(1.01e5,1.2,0,1.21,1.4,-1) = (1.630779226806516, 155686.45, -109.71829086299334)
"""
function after_shock(p1,ρ1,u1,Ms,γ,sdir)
    @assert Ms >= 1
    c1 = sound_speed(rho = ρ1, p = p1, gamma = γ)
    us = Ms * c1 * sdir
    U1 = u1 - us
    M1 = U1 / c1
    r = (γ+1)*M1^2 / (2 + (γ-1)*M1^2)
    U2 = U1 / r
    u2 = U2 + us
    ρ2 = ρ1 * r
    p2 = p1 * (1 + 2*γ/(γ+1) * (M1^2-1))
    # c2 = sp(γ,p2,ρ2)
    # M2 = abs(U2 / c2)
    # M2_predicted = sqrt( (1+(γ-1)/2*M1^2) / (γ*M1^2-(γ-1)/2) )
    # println(M2)
    # println(M2_predicted)
    return ρ2, p2, u2
end

function fill_fluid!(f::Fluid, cell::Cell)
    f.cells = distribute(reshape([copy!(cell) for i in 1:MK.product(f.nmesh .+ f.ng*2)], Tuple(f.nmesh .+ f.ng*2)...), procs = workers(), dist = f.dist)
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.cells)
            bias = [inds[k][1] - 1 for k = 1:f.dim]
            local_inds = Tuple([inds[k] .- bias[k] for k in 1:f.dim])
            local_cells = getindex(localpart(f.cells), local_inds...)
            for i in CartesianIndices(local_cells)
                c = local_cells[i]
                c.i = [i[k] + bias[k] for k = 1:f.dim]
                c.x = f.d .* (c.i .- f.ng .- 0.5) .+ f.point1
            end
        end
    end
    f.background_is_filled = true  
end

function fill_fluid!(f::Fluid, cell::Cell, point1::Array, point2::Array)
    if !f.background_is_filled
        error("Please fill the background before this calling.")
    else
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(f.cells)
                bias = [inds[k][1] - 1 for k = 1:f.dim]
                local_inds = Tuple([inds[k] .- bias[k] for k in 1:f.dim])
                local_cells = getindex(localpart(f.cells), local_inds...)
                for c in local_cells
                    if MK.between(c.x, point1, point2)
                        copy_to!(cell, c)
                    end    
                end
            end
        end        
    end
    
end

