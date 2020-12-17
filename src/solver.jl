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

"""
Don't use deepcopy for DArray. 
"""
function backup_w!(f)
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                localpart(f.wb)[i-bias[1], j-bias[2], k-bias[3]] = f.w[i,j,k]
            end
        end
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
    @sync for pid in workers()
        @spawnat pid begin

            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]

            for i in inds[1], j in inds[2], k in inds[3]
                if MK.betweeneq([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
                    w = coeff[1] * f.w[i,j,k] + coeff[2] * f.wb[i,j,k] + coeff[3] * f.rhs[i,j,k] * dt

                    # w = correct_cell_w(w, f.para["gamma"])

                    localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w

                    rho, u, e, p = w_to_status(w, f.para["gamma"])  

                    localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                    localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                    localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                    localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                end
            end
        end
    end
end

function update_rhs!(f::Fluid)

    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]

            for i in inds[1], j in inds[2], k in inds[3]

                if MK.betweeneq([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
                    
                    rhs = zeros(Float64, 5)
                    ws = Array{Float64,2}(undef, 5, 5)
                    # axis 1
                    for jj in 1:5
                        sw = f.w[i+jj-3, j, k]
                        for ii = 1:5
                            ws[ii, jj] = sw[ii]
                        end
                    end
                    rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 1) / f.d[1]

                    # if (i,j,k) == (17,3,1)
                    #     println("rhs = ", rhs, "  -- x")
                    # end
                    
                    # if rhs[1] != 0.0
                    # println("-- update_rhs: 1 --")
                    # println(rhs)
                    # println(ws)
                    # end

                    # axis 2
                    if f.realdim > 1
                        for jj in 1:5
                            sw = f.w[i, j+jj-3, k]
                            for ii = 1:5
                                ws[ii, jj] = sw[ii]
                            end
                        end
                        rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 2) / f.d[2]  

                        # if (i,j,k) == (17,3,1)
                        #     println("rhs = ", rhs, "  -- y")
                        # end

                        # axis 3
                        if f.realdim > 2
                            for jj in 1:5
                                sw = f.w[i, j, k+jj-3]
                                for ii = 1:5
                                    ws[ii, jj] = sw[ii]
                                end
                            end
                            rhs += get_dflux!(f.para["reconst scheme"], ws, f.para, 3) / f.d[3] 

                            # if (i,j,k) == (17,3,1)
                            #     println("rhs = ", rhs, "  -- z")
                            # end
                        end  
                    end 
                    
                    localpart(f.rhs)[i-bias[1], j-bias[2], k-bias[3]] = rhs

                    # if abs(rhs[1]) > 1e10 && (i,j,k) == (17,3,1)
                    #     println("rhs = ", rhs)
                    #     println((i,j,k))
                    #     exit()
                    # end
                end
            end
        end
    end
end 

function correct_cell_w(w, gamma)
    if sign.([w[1], w[end], pressure(w, gamma)]) == [-1.0, -1.0, -1.0]
        # @warn "Negative mass density"
        # println(c)
    elseif w[1] < 0 || w[end] < 0 || pressure(w, gamma) < 0
        # println("-- correct_cell_w: 1 --")
        # println("c.i = ", c.i)

        w = zeros(Float64, 5)
    end
    return w
end