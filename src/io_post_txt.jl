function save_mesh(f::Fluid, filepath)
    for k in 1:f.dim
        open(filepath*"fluid_x"*string(k)*".txt","w") do file
            writedlm(file, [f.d[k]*(i-0.5-f.ng)+f.point1[k]  for  i=1:f.nmesh[k]+f.ng*2])
        end
    end    
end

function save_to_txt(f::Fluid, field, fname; axis = 1)
    d = getfield(f.cells[1], field)
    if length(d) == 1
        data = Array{typeof(d)}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
        for k in eachindex(f.cells)
            data[k] = getfield(f.cells[k], field)
        end      
    else
        data = Array{typeof(d[1])}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
        for k in eachindex(f.cells)
            data[k] = getfield(f.cells[k], field)[axis]
        end
    end
    open(fname, "w") do file
        writedlm(file, data)
    end
end

function save_to_fig(f::Fluid; dataname::String = "rho", axis::Int = 1, frame::Int = 1, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,1,1])
    if plotdim in ("2D", "2D-x", "2D-y")
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        y = [f.d[2]*(i-0.5-f.ng)+f.point1[2]  for  i=1:f.nmesh[2]+f.ng*2]
        if dataname in ("rho","p","e")
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k], Symbol(dataname))
            end
        elseif dataname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end  
        else
            error("undef dataname")          
        end  
        ioff()
        if plotdim == "2D-x"
            output_plot_func(x, data[:, floor(Int, size(data, 2)/2)], levels=levels)
        elseif plotdim == "2D-y"
            output_plot_func(y, data[floor(Int, size(data, 1)/2), :], levels=levels)  
        else        
            output_plot_func(x, y, data, levels =levels)
        end
    elseif plotdim == "1D"
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        data = Array{Float64}(undef, length(x))
        for k in eachindex(f.cells)
            data[k] = getfield(f.cells[k], Symbol(dataname))
        end
        ioff()
        output_plot_func(x, data, levels=levels)        
    else
        error("undef plotdim")
    end
    savefig(figpath*dataname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end