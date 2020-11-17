function advance!(f::Fluid, dt::Float64)
    
    # println("-- advance_fluid: 1 --")
    # showfield!(f.cells, "rho", 13:20)

    backup_w!(f)

    # println("-- advance_fluid: 2 --")
    # showfield!(f.cells, "rho", 13:20)

    # @warn("rk is limited")
    for rk = 1:3
        # println("-- rk = ",rk, " --")

        # println("-- advance_fluid: 3 --")
        # showfield!(f.cells, "flux", 13:20)

        update_rhs!(f)

        # println("-- advance_fluid: 4 --")
        # showfield!(f.cells, "rhs", 13:20)

        update_cells!(f, rk, dt)

        # println("-- advance_fluid: 5 --")
        # showfield!(f.cells, "rho", 13:20)

        # 下面这种写法是错的，因为wb不会随着rk更新。
        # update_fluxes_and_cells!(f, rk = rk, dt = dt)
        update_boundaries!(f)

        # println("-- advance_fluid: 6 --")
        # showfield!(f.cells, "rho", 13:20)

        
    end

end

function solve!(f::Fluid; CFL::Float64 = 0.3, maxtime::Float64 = 1, maxframe::Int = 1, cutframe::Int = 1, varname::String = "rho", axis::Int = 1, filepath::String = "../outputdata/", draw::Bool = false, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,0,0])
    time = 0.
    frame = 0
    dt = 0.
    println("|||||||||||||||||||||||||||||||||||||||||||")
    @printf "Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100

    if OUTPUTDATA
        output!(Int(frame/cutframe), time, filepath = filepath)
        output!(f, varname = varname, axis = axis, frame = Int(frame/cutframe), filepath = filepath)
        if draw
            outputfig!(f, frame = Int(frame/cutframe), varname = varname, axis = axis, figpath = figpath, plotdim = plotdim, levels=levels)
        end
    end
    while time < maxtime && frame < maxframe
        dt = time_step!(f, CFL = CFL)
        advance!(f, dt)
        time += dt
        frame += 1
        print("Current frame = ",frame,"\r")
        if frame%cutframe == 0
            println("|||||||||||||||||||||||||||||||||||||||||||")
            @printf "Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100
            if OUTPUTDATA
                output!(Int(frame/cutframe), time, filepath = filepath)
                output!(f, varname = varname, axis = axis, frame = Int(frame/cutframe), filepath = filepath)
                if draw
                    outputfig!(f, frame = Int(frame/cutframe), varname = varname, axis = axis, figpath = figpath, plotdim = plotdim,levels=levels)
                end                
            end
        end
    end
end

"""
像time_step!()这样仅读取的函数，则无需远程，可直接读取DArray。
"""
function time_step!(f::Fluid; CFL::Float64 = 0.1)
    smax = 0.
    @sync for pid in workers()
        smaxtmptmp= @fetchfrom pid begin
            smaxtmp = 0.
            for c in localpart(f.cells)
                if c.rho == 0.
                    s = 0.
                else
                    s = sound_speed(rho = c.rho, p = c.p, gamma = f.constants["gamma"])
                end
                
                smaxtmp = maximum([smaxtmp; s .+ map(abs, c.u)])
            end
            smaxtmp
        end
        smax = max(smax, smaxtmptmp)
    end
    dt = minimum(f.d) / smax * CFL
    if dt < TOL_STEP
        error( "Too small time step!")
    end
    return dt
end

function update_cells!(f::Fluid, rk::Int, dt::Float64)
    coeff = RK_COEFF[:, rk]
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                if MK.betweeneq(c.x, f.point1, f.point2)
            
                    c.w = coeff[1] * c.w + coeff[2] * c.wb + coeff[3] * c.rhs * dt


                    # correct_cell_w!(c, f.constants["gamma"])
                    w2states!(c, f.constants["gamma"])  

                end
            end
        end
    end
end

function update_rhs!(f::Fluid)
    if f.dim == 1
        @sync for pid in workers()
            @spawnat pid begin
                ws = Array{Float64,2}(undef, 2+f.dim, 5)
    
                for c in localpart(f.cells) # 千万别漏写localpart 
                    if MK.betweeneq(c.x, f.point1, f.point2)
                        c.rhs = zeros(Float64, 2+f.dim)
                        ci = copy(c.i)
    
                        for j in 1:5
                            sw = f.cells[ci[1]+j-3].w
                            for i = 1:2+f.dim
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 1) / f.d[1] 
                                         
                    end
                end
            end
        end
    elseif f.dim == 2
        @sync for pid in workers()
            @spawnat pid begin
                ws = Array{Float64,2}(undef, 2+f.dim, 5)
    
                for c in localpart(f.cells) # 千万别漏写localpart 
                    if MK.betweeneq(c.x, f.point1, f.point2)
                        c.rhs = zeros(Float64, 2+f.dim)
                        ci = copy(c.i)
    
                        for j in 1:5
                            sw = f.cells[ci[1]+j-3,ci[2]].w
                            for i = 1:4
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 1) / f.d[1] 
                        
    
                        for j in 1:5
                            sw = f.cells[ci[1],ci[2]+j-3].w
                            for i = 1:4
                                ws[i,j] = sw[i]
                            end
                        end    
    
                        c.rhs += get_dflux!(f.reconst_scheme, ws, f.constants, 2) / f.d[2]                    
                    end
                end
            end
        end
    else
        error("undef dim")
    end
end 

# """
# 这样写会慢30%
# """
# function testupdate_rhs!(f::Fluid)
#     @sync for pid in workers()
#         @spawnat pid begin
#             ws = Array{Float64,2}(undef, 2+f.dim, 5)
#             for c in localpart(f.cells) # 千万别漏写localpart 
#                 if MK.betweeneq(c.x, f.point1, f.point2)     
#                     c.rhs = zeros(Float64, 2+f.dim)
#                     for axis = 1:f.dim
#                         i = copy(c.i)
                         
#                         for k = 1:5
#                             i[axis] = c.i[axis] + k - 3
#                             ws[:,k] = f.cells[Tuple(i)...].w
#                         end 
                        
#                         c.rhs += (get_flux!(f.reconst_scheme, ws[:,1:4], f.constants, axis = axis) - get_flux!(f.reconst_scheme, ws[:,2:5], f.constants, axis = axis)) / f.d[axis]
#                     end                       
                    
#                 end
#             end
#         end
#     end
# end 

