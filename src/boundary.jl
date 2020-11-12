
"""
gcells是要被修改的边界外单元。cells是作为参照的边界内单元。
"""
function ghost_boundary_cell!(cells::SubArray, gcells::Array{Cell}, boundary::String; axis::Int = -10)
    copy_to!(cells, gcells)
    if boundary == "refl"
        factor = - 1
    else
        factor = 1
    end

    for k in eachindex(cells)

        gcells[k].u[axis] *= factor
        gcells[k].w[axis + 1] *= factor
    end 
end

function set_bounds!(f::Fluid, boundaries::Array{String})
    f.boundaries = boundaries
    update_boundaries!(f)
end


"""
一个合理的思路是：
1. 把要修改的子数组的索引modrange发送到pid分区。

2. 在pid分区上操作。首先确定要修改的子数组索引modrange与pid分区索引pidrange的交集intrange。

3. 对于不需要修改的子数组，直接从本地拉取到pid分区。

4. 其余不需要修改的参数也从本地拉取即可。

5. 在pid分区调用函数计算，最后把结果拉取到本地。
"""
function update_boundaries!(f::Fluid)
    nsize = (f.nmesh .+ f.ng*2) ./f.dist
    for n in nsize
        if n < f.ng * 2
            error("Too small nsize")
        end
    end
    @sync for pid in workers() 
        @spawnat pid begin
            pid_range = localindices(f.cells)
            bias = [r[1] for r in pid_range] .- 1
            for axis in 1:f.dim
                for side in 1:2
                    for i in (f.ng+f.nmesh[axis])*(side-1) .+ (1:f.ng)
                        int_range = map(intersect, pid_range, Tuple([k == axis ? (i:i) : (1:size(f.cells, k)) for k in 1:f.dim]))
                        if length_of_tuples(int_range) > 0  # 如果intrange索引的“覆盖面积”是0,那么就不做计算。
                            ref_i = 2*(f.ng+f.nmesh[axis]*(side-1))+1-i
                            ref_range = Tuple([k == axis ? (ref_i:ref_i) : int_range[k] for k = 1:f.dim])
                            local_int_range = Tuple([int_range[k] .- bias[k] for k in 1:f.dim])
                            ghost_boundary_cell!(f.cells[ref_range...], localpart(f.cells)[local_int_range...], f.boundaries[axis, side], axis = axis)   
                            # println(("side",side,"pid_range = ",pid_range,"bias = ",bias,"int_range = ",int_range,"local_int_range",local_int_range,"ref_range = ",ref_range))                                                
                        end
                    end
                end
            end
        end
    end
    # showfield!(f.cells, "rho")
end