function copy_fluid!(f::Fluid)
    f1 = Fluid(realdim = f.realdim, point1 = f.point1, point2 = f.point2, nmesh = f.nmesh, ng = f.ng, para = f.para)

    for id in CartesianIndices(f.rho)
        f1.rho[id] = f.rho[id]
        f1.e[id] = f.e[id]
        f1.p[id] = f.p[id]
        f1.mark[id] = f.mark[id]

        f1.u[:,id] = f.u[:,id]
        f1.w[:,id] = f.w[:,id]
        f1.wb[:,id] = f.wb[:,id]
        f1.rhs[:,id] = f.rhs[:,id]
        f1.target_id[:,id] = f.target_id[:,id]
    end

    f1.boundx = f.boundx
    f1.boundy = f.boundy
    f1.boundz = f.boundz

    return f1
end