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

        for filled_values in f.para["fill forever"]
            fill_fluid!(f, filled_values...)
        end

        update_bounds!(f)

        # println("-- advance_fluid: 6 --")
        # showfield!(f.cells, "rho", 13:20)

        
    end
end

function backup_w!(f)
    @sync @distributed for id in eachindex(f.w)
        f.wb[id] = f.w[id]
    end     
end

function time_step!(f::Fluid; CFL::Float64 = 0.1)
    
    # correct_cells!(f)

    smax = 0.
    smax = @sync @distributed (max) for id in eachindex(f.rho)
        s = sound_speed(f.rho[id], f.p[id], f.para["gamma"])
        u = f.u[id]
        smax = maximum([smax; s .+ map(abs, u)])
    end
    dt = minimum(f.d[1:f.realdim]) / smax * CFL
    if isnan(dt)
        # if smax > 0.
            println("smax = ", smax)
            error("NaN")
        # end
    end
    # if dt < TOL_STEP
    #     error( "Too small time step!")
    # end
    return dt
end

function update_cells!(f::Fluid, rk::Int, dt::Float64)
    coeff = RK_COEFF[:, rk]

    @sync @distributed for id in CartesianIndices(f.rho)
        if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
        # if f.marker[id] > 0 && MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
                
                w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * f.rhs[:, id] * dt

                f.w[:, id] = w
    
                rho, u, e, p = w_to_status(w, f.para["gamma"]) 
    
                f.rho[id] = rho
                f.u[:,id] = u
                f.e[id] = e
                f.p[id] = p
            end
    end

    correct_cells!(f)
end

function correct_cells!(f::Fluid)
    @sync @distributed for id in CartesianIndices(f.rho)
        if MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
        # if f.marker[id] > 0 && MK.betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)

            w = correct_cell_w(f.w[:, id], f.para["gamma"], f.para["rho0"], f.para["u0"], f.para["e0"])  ## This correction may cause mass loss.
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
    @sync @distributed for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if MK.betweeneq([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
        # if f.marker[id] > 0 && MK.betweeneq([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
                rhs = zeros(Float64, 5)
                ws = Array{Float64,2}(undef, 5, 5)
                # axis 1
                for jj in 1:5
                    ws[:, jj] = f.w[:, i+jj-3, j, k]
                end
                rhs += get_dflux!(id, f.para["reconst scheme"], ws, f.para, 1) / f.d[1]
    
                # axis 2
                if f.realdim > 1
                    for jj in 1:5
                        ws[:, jj] = f.w[:, i, j+jj-3, k]
                    end
                    rhs += get_dflux!(id, f.para["reconst scheme"], ws, f.para, 2) / f.d[2]  
    
                    # axis 3
                    if f.realdim > 2
                        for jj in 1:5
                            ws[:, jj] = f.w[:, i, j, k+jj-3]
                        end
                        rhs += get_dflux!(id, f.para["reconst scheme"], ws, f.para, 3) / f.d[3] 
                    end  
                end 
                
                f.rhs[:,id] = rhs
                   
        end
    end    
end 