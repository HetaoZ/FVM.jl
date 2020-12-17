function advance!(f::Fluid, dt::Float64)
    
    # println("-- advance_fluid: 1 --")
    # showfield!(f.cells, "rho", 13:20)

    backup_w!(f)

    # println("-- advance_fluid: 2 --")
    # showfield!(f.cells, "rho", 13:20)

    # @warn("rk is limited")
    for rk = 1:3

        # println("-- advance_fluid: 3 --")
        # showfield!(f.cells, "flux", 13:20)

        update_rhs!(f)

        # println("-- advance_fluid: 4 --")
        # showfield!(f.cells, "rhs", 13:20)

        update_cells!(f, rk, dt)

        # println("-- advance_fluid: 5 --")
        # showfield!(f.cells, "rho", 13:20)

        update_bounds!(f)

        # println("-- advance_fluid: 6 --")
        # showfield!(f.cells, "rho", 13:20)

        
    end
end

function backup_w!(f)
    for id in eachindex(f.w)
        f.wb[id] = f.w[id]
    end     
end

function time_step!(f::Fluid; CFL::Float64 = 0.1)
    smax = 0.
    for id in eachindex(f.rho)
        s = sound_speed(f.rho[id], f.p[id], f.para["gamma"])
        u = f.u[id]
        smax = maximum([smax; s .+ map(abs, u)])
        
    end
    dt = minimum(f.d[1:f.realdim]) / smax * CFL
    if isnan(dt)
        # if smax > 0.
            println("smax = ", smax)
            exit()
        # end
    end
    # if dt < TOL_STEP
    #     error( "Too small time step!")
    # end
    return dt
end

function update_cells!(f::Fluid, rk::Int, dt::Float64)
    coeff = RK_COEFF[:, rk]

    for id in CartesianIndices(f.rho)
        if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
            w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * f.rhs[:, id] * dt

            ## This correction may cause mass loss.
            # w = correct_cell_w(w, f.para["gamma"]) 

            f.w[:, id] = w

            rho, u, e, p = w_to_status(w, f.para["gamma"]) 

            f.rho[id] = rho
            f.u[:,id] = u
            f.e[id] = e
            f.p[id] = p
        end
    end
end

function update_rhs!(f::Fluid)
    for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if MK.betweeneq([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
            rhs = zeros(Float64, 5)
            ws = Array{Float64,2}(undef, 5, 5)
            # axis 1
            for jj in 1:5
                sw = f.w[:, i+jj-3, j, k]
                for ii = 1:5
                    ws[ii, jj] = sw[ii]
                end
            end
            rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 1) / f.d[1]

            # axis 2
            if f.realdim > 1
                for jj in 1:5
                    sw = f.w[:, i, j+jj-3, k]
                    for ii = 1:5
                        ws[ii, jj] = sw[ii]
                    end
                end
                rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 2) / f.d[2]  

                # axis 3
                if f.realdim > 2
                    for jj in 1:5
                        sw = f.w[:, i, j, k+jj-3]
                        for ii = 1:5
                            ws[ii, jj] = sw[ii]
                        end
                    end
                    rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 3) / f.d[3] 
                end  
            end 
            
            f.rhs[:,id] = rhs
        end
    end    
end 

function correct_cell_w(w, gamma)
    if sign.([w[1], w[end], pressure(w, gamma)]) == [-1.0, -1.0, -1.0]
        # @warn "Negative mass density"
        # println(w)
    elseif w[1] < 0 || w[end] < 0 || pressure(w, gamma) < 0
        w = zeros(Float64, 5)
    end
    return w
end