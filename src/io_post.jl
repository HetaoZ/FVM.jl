function meshxy(x::Array, y::Array; dim::Int = 1)
    m=maximum(size(x))
    n=maximum(size(y))
    X=Array{Float64}(undef,m,n)
    Y=Array{Float64}(undef,m,n)
    for i=1:m,j=1:n
        X[i,j]=x[i]
        Y[i,j]=y[j]
    end
    if dim == 1
        return reshape(X,m*n,1),reshape(Y,m*n,1)
    elseif dim == 2
        return X,Y
    else
        error( "The dimension must be <= 2")
    end
end

function read_frame(n::Int; varname::String = "rho", filepath::String = "../outputdata/")
    var = readdlm(filepath*"fluid_"*varname*"_"*string(n+FRAME_BASE)*".txt")
    return var
end

function read_mesh(dim::Int, n::Int; varname::String = "rho", filepath::String = "../outputdata/")
    return Tuple([readdlm(filepath*"fluid_x"*string(k)*".txt") for k in 1:dim])
end

function output!(frame::Int, time::Float64; filepath::String = "../outputdata/")
    if frame == 0
        open(filepath*"time.txt","w") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    else
        open(filepath*"time.txt","a") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    end
end

"""
axis表示输出变量的方向。
"""
function output!(f::Fluid; varname::String = "rho", axis::Int = 1, frame::Int = 1, filepath::String = "../outputdata/")
    if varname == "mesh"
        for k in 1:f.dim
            open(filepath*"fluid_x"*string(k)*".txt","w") do file
                writedlm(file, [f.d[k]*(i-0.5-f.ng)+f.point1[k]  for  i=1:f.nmesh[k]+f.ng*2])
            end
        end        

    else
        # essential vars
        for essential_varname in ("rho", "p")
            data = Array{typeof(getfield(f.cells[1],Symbol(essential_varname)))}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k],Symbol(essential_varname))
            end      
            open(filepath*"fluid_"*essential_varname*"_"*string(frame+FRAME_BASE)*".txt","w") do file
                writedlm(file, data)
            end           
        end           
        # reviewed vars
        if varname in ("e",)
            data = Array{typeof(getfield(f.cells[1],Symbol(varname)))}(undef, Tuple(f.nmesh .+ 2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k],Symbol(varname))
            end
        end
        if varname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end      
            open(filepath*"fluid_"*varname*"_"*string(frame+FRAME_BASE)*".txt","w") do file
                writedlm(file, data)
            end      
        end

    end
end

function output_plot_func(x::Array{Float64}, data::Array{Float64}; levels::Vector = [0, 1, 1])
    figure(1, figsize=(8, 6))
    c = plot(x, data)

    # ref
    A = readdlm("../../../CommonData/monasse2012-piston-p.txt")
    plot(A[:,1], A[:,2], color = "r")
    
    ylim(levels[1], levels[2])
end

function output_plot_func(x::Array{Float64}, y::Array{Float64}, data::Array{Float64}; levels::Vector = [0, 1, 1], cmap::String = "jet", cbar::Bool = true)
    scale = FIGHEIGHT / (y[end] - y[1])
    figure(1, figsize=((x[end] - x[1]) * scale + 4 * cbar, FIGHEIGHT))
    X, Y = meshxy(x, y, dim = 2)
    if levels[2] <= levels[1] || levels[3] < 2
        levels = (minimum(data), maximum(data), 30)
    end
    c = contourf(X, Y, data, levels=LinRange(Tuple(levels)...), cmap=cmap);
    c1 = contour(X, Y, data, levels=LinRange(Tuple(levels)...));
    # c2 = scatter(meshxy(x, y, dim = 1)..., marker = "+");
    xlim([x[3], x[end-2]]); ylim([y[3], y[end-2]])
    if cbar 
        colorbar(c)
    end
end

function outputfig(frame::Int; varname::String = "rho", filepath::String = "outputdata/", figpath::String = "outputfig/",  levels::Vector = [0, 1, 1], cmap::String = "jet", cbar::Bool = true, plotdim::String = "2D")
    if plotdim == "1D"
        x = read_mesh(1, frame, varname = varname, filepath = filepath)
        data = read_frame(frame, varname = varname, filepath = filepath)
        ioff()
        output_plot_func(x[1], data, levels=levels)
    elseif plotdim in ("2D", "2D-x", "2D-y")
        x, y = read_mesh(2, frame, varname = varname, filepath = filepath)
        data = read_frame(frame, varname = varname, filepath = filepath)
        ioff()
        if plotdim == "2D-x"
            output_plot_func(x, data[:, floor(Int, size(data, 2)/2)], levels=levels)
        elseif plotdim == "2D-y"
            output_plot_func(y, data[floor(Int, size(data, 1)/2), :], levels=levels)  
        else
            output_plot_func(x, y, data, levels = levels, cmap = cmap, cbar = cbar)    
        end 
    else
        error("undef dim")       
    end
    savefig(figpath*varname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end

function outputfig(frames::AbstractRange; varname::String = "rho", filepath::String = "outputdata/", figpath::String = "outputfig/",  levels::Array = [0, 0, 0], cmap::String = "jet", cbar::Bool = true, plotdim::String = "2D")
    for frame in frames
        outputfig(frame, varname = varname, filepath = filepath, figpath = figpath, levels = levels, cmap = cmap, cbar = cbar, plotdim = plotdim)
    end
end

"""
从内存中读取远比硬盘快。但也需要写一个从硬盘读取并绘图的函数。
"""
function outputfig!(f::Fluid; varname::String = "rho", axis::Int = 1, frame::Int = 1, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,1,1])
    if plotdim in ("2D", "2D-x", "2D-y")
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        y = [f.d[2]*(i-0.5-f.ng)+f.point1[2]  for  i=1:f.nmesh[2]+f.ng*2]
        if varname in ("rho","p","e")
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k], Symbol(varname))
            end
        elseif varname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end  
        else
            error("undef varname")          
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
            data[k] = getfield(f.cells[k], Symbol(varname))
        end
        ioff()
        output_plot_func(x, data, levels=levels)        
    else
        error("undef plotdim")
    end
    savefig(figpath*varname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end