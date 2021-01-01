
"""
Corners are not to be updated
"""
function update_bounds!(f::Fluid)
    ng = f.ng
    nmesh = f.nmesh

    @sync @distributed for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if i < ng + 1
            marker = f.marker[2*ng+1-i, j, k]
            rho = f.rho[2*ng+1-i, j, k]
            u = f.u[:, 2*ng+1-i, j, k]
            u[1] *= BOUND_TYPE[f.boundx[1]]
            e = f.e[2*ng+1-i, j, k]
            p = f.p[2*ng+1-i, j, k]
            w = status_to_w(rho, u, e)

            f.marker[id] = marker
            f.rho[id] = rho
            f.u[:, id] = u
            f.e[id] = e
            f.p[id] = p
            f.w[:, id] = w
        end
        if i > nmesh[1] + ng
            marker = f.marker[2*(nmesh[1]+ng)+1-i, j, k]
            rho = f.rho[2*(nmesh[1]+ng)+1-i, j, k]
            u = f.u[:, 2*(nmesh[1]+ng)+1-i, j, k]
            u[1] *= BOUND_TYPE[f.boundx[2]]
            e = f.e[2*(nmesh[1]+ng)+1-i, j, k]
            p = f.p[2*(nmesh[1]+ng)+1-i, j, k]
            w = status_to_w(rho, u, e)

            f.marker[id] = marker
            f.rho[id] = rho
            f.u[:, id] = u
            f.e[id] = e
            f.p[id] = p
            f.w[:, id] = w
        end
        
        if f.realdim > 1
            if j < ng + 1
                marker = f.marker[i, 2*ng+1-j, k]
                rho = f.rho[i, 2*ng+1-j, k]
                u = f.u[:, i, 2*ng+1-j, k]
                u[2] *= BOUND_TYPE[f.boundy[1]]
                e = f.e[i, 2*ng+1-j, k]
                p = f.p[i, 2*ng+1-j, k]
                w = status_to_w(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
            if j > nmesh[2] + ng
                marker = f.marker[i, 2*(nmesh[2]+ng)+1-j, k]
                rho = f.rho[i, 2*(nmesh[2]+ng)+1-j, k]
                u = f.u[:, i, 2*(nmesh[2]+ng)+1-j, k]
                u[2] *= BOUND_TYPE[f.boundy[2]]
                e = f.e[i, 2*(nmesh[2]+ng)+1-j, k]
                p = f.p[i, 2*(nmesh[2]+ng)+1-j, k]
                w = status_to_w(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w                         
            end
            if f.realdim > 2
                if k < ng + 1
                    marker = f.marker[i, j, 2*ng+1-k]
                    rho = f.rho[i, j, 2*ng+1-k]
                    u = f.u[:, i, j, 2*ng+1-k]
                    u[3] *= BOUND_TYPE[f.boundz[1]]
                    e = f.e[i, j, 2*ng+1-k]
                    p = f.p[i, j, 2*ng+1-k]
                    w = status_to_w(rho, u, e)

                    f.marker[id] = marker
                    f.rho[id] = rho
                    f.u[:, id] = u
                    f.e[id] = e
                    f.p[id] = p
                    f.w[:, id] = w                     
                end
                if k > nmesh[3] + ng
                    marker = f.marker[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                    rho = f.rho[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                    u = f.u[:, i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                    u[3] *= BOUND_TYPE[f.boundz[2]]
                    e = f.e[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                    p = f.p[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                    w = status_to_w(rho, u, e)

                    f.marker[id] = marker
                    f.rho[id] = rho
                    f.u[:, id] = u
                    f.e[id] = e
                    f.p[id] = p
                    f.w[:, id] = w                     
                end
            end
        end
    end 
end