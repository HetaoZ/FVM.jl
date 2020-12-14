
"""
Corners are not to be updated.
"""
function update_bounds!(f::Fluid)
    ng = f.ng
    nmesh = f.nmesh
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            # println("inds = ", inds)
            # println("bias = ", bias)
            for i in inds[1], j in inds[2], k in inds[3]
                if i < f.ng+1
                    rho = f.rho[2*ng+1-i, j, k]
                    u = f.u[2*ng+1-i, j, k]
                    u[1] *= BOUND_TYPE[f.boundx[1]]
                    e = f.e[2*ng+1-i, j, k]
                    p = f.p[2*ng+1-i, j, k]
                    w = status_to_w(rho, u, e)
                    
                    localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                    localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                    localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                    
                    localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                    
                    localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w
                end                        
                if i > nmesh[1]+ng
                    rho = f.rho[2*(nmesh[1]+ng)+1-i, j, k]
                    u = f.u[2*(nmesh[1]+ng)+1-i, j, k]
                    u[1] *= BOUND_TYPE[f.boundx[2]]
                    e = f.e[2*(nmesh[1]+ng)+1-i, j, k]
                    p = f.p[2*(nmesh[1]+ng)+1-i, j, k]
                    w = status_to_w(rho, u, e)

                    localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                    localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                    localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                    
                    localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                    
                    localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w                        
                end 
                if f.realdim > 1
                    if j < f.ng+1
                        rho = f.rho[i, 2*ng+1-j, k]
                        u = f.u[i, 2*ng+1-j, k]
                        u[2] *= BOUND_TYPE[f.boundy[1]]
                        e = f.e[i, 2*ng+1-j, k]
                        p = f.p[i, 2*ng+1-j, k]
                        w = status_to_w(rho, u, e)

                        localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                        localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                        localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                        
                        localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                        
                        localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w       
                    end
                    if j > nmesh[2]+ng
                        rho = f.rho[i, 2*(f.nmesh[2]+f.ng)+1-j, k]
                        u = f.u[i, 2*(f.nmesh[2]+f.ng)+1-j, k]
                        u[2] *= BOUND_TYPE[f.boundy[2]]
                        e = f.e[i, 2*(f.nmesh[2]+f.ng)+1-j, k]
                        p = f.p[i, 2*(f.nmesh[2]+f.ng)+1-j, k]
                        w = status_to_w(rho, u, e)

                        localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                        localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                        localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                        
                        localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                        
                        localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w                            
                    end
                    if f.realdim > 2
                        if k < ng+1 || k > nmesh[3]+ng
                            if k < f.ng+1
                                rho = f.rho[i, j, 2*ng+1-k]
                                u = f.u[i, j, 2*ng+1-k]
                                u[3] *= BOUND_TYPE[f.boundz[1]]
                                e = f.e[i, j, 2*ng+1-k]
                                p = f.p[i, j, 2*ng+1-k]
                                w = status_to_w(rho, u, e)
                            else
                                rho = f.rho[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                                u = f.u[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                                u[3] *= BOUND_TYPE[f.boundz[2]]
                                e = f.e[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                                p = f.p[i, j, 2*(f.nmesh[3]+f.ng)+1-k]
                                w = status_to_w(rho, u, e)
                            end
                            localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = rho
                            localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = u
                            localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = p
                            
                            localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = e
                            
                            localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = w
                        end
                    end
                end

                              
            end
        end
    end   
end