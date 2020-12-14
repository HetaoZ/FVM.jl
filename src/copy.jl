function copy_fluid!(f::Fluid)
    f1 = Fluid(realdim = f.realdim)

    f1.point1 = copy(f.point1)
    f1.point2 = copy(f.point2)
    f1.nmesh = copy(f.nmesh)
    f1.d = copy(f.d)
    f1.ng = f.ng

    f1.dist = copy(f.dist)
    f1.ndiv = f.ndiv

    f1.NX = f.NX
    f1.NY = f.NY
    f1.NZ = f.NZ

    f1.x = f.x
    f1.y = f.y
    f1.z = f.z

    f1.rho = copy_darray!(f.rho, f.ndiv, f.dist)
    f1.u = copy_darray!(f.u, f.ndiv, f.dist)
    f1.e = copy_darray!(f.e, f.ndiv, f.dist)
    f1.p = copy_darray!(f.p, f.ndiv, f.dist)

    f1.w = copy_darray!(f.w, f.ndiv, f.dist)

    f1.boundx = copy(f.boundx)
    f1.boundy = copy(f.boundy)
    f1.boundz = copy(f.boundz)

    f1.para = f.para
    return f1
end

function copy_darray!(a::DArray{T, N, A} where T where N where A, ndiv::Int, dist::Array{Int}) 
    tmp = Array{eltype(a)}(undef, size(a))
    for i in eachindex(tmp)
        tmp[i] = copy(a[i])
    end
    return distribute(tmp, procs = workers(), dist = dist)
end