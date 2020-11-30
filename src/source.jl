get_sxx(hx::Float64, u::Vector{Float64}, hy::Float64, v::Vector{Float64}, mu::Float64) = 2/3*mu*(2*(u[3] - u[1])/(2*hx) - (v[3] - v[1])/(2*hy))

get_sxy(hy::Float64, u::Vector{Float64}, hx::Float64, v::Vector{Float64}, mu::Float64) = mu*((u[3] - u[1])/(2*hy) + (v[3] - v[1])/(2*hx))

get_syy(hx::Float64, u::Vector{Float64}, hy::Float64, v::Vector{Float64}, mu::Float64) = 2/3*mu*(- (u[3] - u[1])/(2*hx) + 2*(v[3] - v[1])/(2*hy))

function get_vis_item!(u::Array{Float64}, v::Array{Float64}, T::Array{Float64}, para::Dict, d::Vector, rho::Float64)
    if rho == 0.0
        return zeros(Float64, length(d) + 2)
    end
    Re = reynolds_number(rho, para)
    mu = para["mu"]
    beta = - mu/(para["Pr"]*(para["gamma"]-1.0))
    if length(d) == 2
        hx, hy = d[1], d[2]

        # 现在在5x5的点阵里考虑. 编号1~5，中心是3.
        sxx_43 = get_sxx(hx, u[3:5, 3], hy, v[4, 2:4], mu)
        sxx_23 = get_sxx(hx, u[1:3, 3], hy, v[2, 2:4], mu)
    
        sxy_43 = get_sxy(hy, u[4, 2:4], hx, v[3:5, 3], mu) # 注意应力张量对称
        sxy_23 = get_sxy(hy, u[2, 2:4], hx, v[1:3, 3], mu) # 注意应力张量对称
    
        sxy_34 = get_sxy(hy, u[3, 3:5], hx, v[2:4, 4], mu)
        sxy_32 = get_sxy(hy, u[3, 1:3], hx, v[2:4, 2], mu)
    
        syy_34 = get_syy(hx, u[2:4, 4], hy, v[3, 3:5], mu)
        syy_32 = get_syy(hx, u[2:4, 2], hy, v[3, 1:3], mu)
    
        qx_43 = beta * (T[5, 3] - T[3, 3])/(2*hx)
        qx_23 = beta * (T[3, 3] - T[1, 3])/(2*hx)
        
        qy_34 = beta * (T[3, 5] - T[3, 3])/(2*hy)
        qy_32 = beta * (T[3, 3] - T[3, 1])/(2*hy)
    
        bx_43 = u[4, 3]*sxx_43 + v[4, 3]*sxy_43 - qx_43
        bx_23 = u[2, 3]*sxx_23 + v[2, 3]*sxy_23 - qx_23
    
        by_34 = u[3, 4]*sxy_34 + v[3, 4]*syy_34 - qy_34
        by_32 = u[3, 2]*sxy_32 + v[3, 2]*syy_32 - qy_32
        
        vis_item_x = [0.0, sxx_43 - sxx_23, sxy_43 - sxy_23, bx_43 - bx_23] ./ (2*hx) 
        vis_item_y = [0.0, sxy_34 - sxy_32, syy_34 - syy_32, by_34 - by_32] ./ (2*hy) 
    
        vis_item = (vis_item_x .+ vis_item_y) ./ Re
        if any(map(isnan, vis_item))
            println("-- get_vis_item: 0 --")
            display(vis_item_x)
            display(vis_item_y)
            display(Re)
            display([sxx_43 sxx_23 sxy_43 sxy_23 bx_43 bx_23 hx;
                     sxy_34 sxy_32 syy_34 syy_32 by_34 by_32 hy])
            error("get_vis_item error")
        end
    else
        error("get_vis_item error")
    end
    return vis_item
end

