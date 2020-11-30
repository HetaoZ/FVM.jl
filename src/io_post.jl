function review(f::Fluid)
    println("-- short summary of Fluid --")
    println("# parameters")
    for k in keys(f.para)
        println("  ",k," : ",f.para[k])
    end
    println("# mesh")
    println("  dimension : ", f.dim)
    println("  domain : ", Tuple(f.point1)," -> ", Tuple(f.point2))
    println("  number of cells : ", length(f.cells))
    println("  size of cells : ", size(f.cells))
    println("# distribution")
    showdist!(f.cells)
    println("# boundaries")
    for k = 1:f.dim
        println("  axis "*string(k)*" : ", (f.boundaries[k,1], f.boundaries[k,2]))
    end
    println("# states")
    rho = fetch_data(f, :rho)
    p = fetch_data(f, :p)
    u = fetch_data(f, :u)
    println("  rho : ", [minimum(rho), maximum(rho)])
    println("  p : ", [minimum(p), maximum(p)])
    println("  u : ", [minimum(u), maximum(u)])
    println()
end

function save_review(f, fname)
    open(fname*".txt", "w") do file
        write(file,"-- short summary of Fluid --")
        write(file,"# parameters")
        for k in keys(f.para)
            write(file,"  "*string(k)*" : "*string(f.para[k]))
        end
        write(file,"# mesh")
        write(file,"  dimension : "*string(f.dim))
        write(file,"  domain : "*string(Tuple(f.point1))*" -> "*string(Tuple(f.point2)))
        write(file,"  number of cells : "*string(length(f.cells)))
        write(file,"  size of cells : "*string(size(f.cells)))
        write(file,"# distribution")
        write(file, nworkers())
        write(file,"# boundaries")
        for k = 1:f.dim
            write(file,"  axis "*string(k)*" : "*string((f.boundaries[k,1], f.boundaries[k,2])))
        end
        write(file,"# states")
        rho = fetch_data(f, :rho)
        p = fetch_data(f, :p)
        u = fetch_data(f, :u)
        write(file,"  rho : "*string([minimum(rho), maximum(rho)]))
        write(file,"  p : "*string([minimum(p), maximum(p)]))
        write(file,"  u : "*string([minimum(u), maximum(u)]))
        write(file,"")
    end
end

function fetch_data(f, field)
    datadim = length(getfield(f.cells[1], field))
    if datadim == 1
        a = Array{typeof(getfield(f.cells[1], field))}(undef, size(f.cells))
        for i in eachindex(f.cells)
            a[i] = getfield(f.cells[i], field)
        end
    else
        a = Array{eltype(getfield(f.cells[1], field))}(undef, datadim, size(f.cells)...)
        for i in eachindex(f.cells)
            a[:, i] = getfield(f.cells[i], field)
        end
    end
    return a
end