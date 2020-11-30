function save_to_vtk(f, datanames, fields, fname)
    vtkfile = create_vtkfile(f, fname)
    for i in eachindex(datanames)
        vtkfile[datanames[i]] = fetch_data(f, fields[i])
    end
    outfiles = vtk_save(vtkfile)
end

function create_vtkfile(f, fname)
    if f.dim == 1
        x = f.point1[1]-(f.ng-0.5)*f.d[1]:f.d[1]:f.point2[1]+(f.ng-0.5)*f.d[1]
        file = vtk_grid(fname, x)
    elseif f.dim == 2
        x = f.point1[1]-(f.ng-0.5)*f.d[1]:f.d[1]:f.point2[1]+(f.ng-0.5)*f.d[1]
        y = f.point1[2]-(f.ng-0.5)*f.d[2]:f.d[2]:f.point2[2]+(f.ng-0.5)*f.d[2]
        file = vtk_grid(fname, x, y)
    else
        x = f.point1[1]-(f.ng-0.5)*f.d[1]:f.d[1]:f.point2[1]+(f.ng-0.5)*f.d[1]
        y = f.point1[2]-(f.ng-0.5)*f.d[2]:f.d[2]:f.point2[2]+(f.ng-0.5)*f.d[2]
        z = f.point1[3]-(f.ng-0.5)*f.d[3]:f.d[3]:f.point2[3]+(f.ng-0.5)*f.d[3]
        file = vtk_grid(fname, x, y, z)
    end
    return file
end

function fetch_data(f, field)
    datadim = length(getfield(f.cells[1], field))
    if datadim == 1
        a = Array{typeof(getfield(f.cells[1], field))}(undef, size(f.cells))
        for i in eachindex(f.cells)
            a[i] = getfield(f.cells[i], field)
        end
    else
        a = Array{eltype(getfield(f.cells[1], field))}(undef, datadim, size(f.cells))
        for i in eachindex(f.cells)
            a[:, i] = getfield(f.cells[i], field)
        end
    end
    return a
end

function save_fluid_mesh(f::Fluid, filepath)
    for k in 1:f.dim
        open(filepath*"fluid_x"*string(k)*".txt","w") do file
            writedlm(file, [f.d[k]*(i-0.5-f.ng)+f.point1[k]  for  i=1:f.nmesh[k]+f.ng*2])
        end
    end    
end