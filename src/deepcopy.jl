"""
deepcopy?
"""
function backup_w!(f::Fluid)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                c.wb = copy(c.w)
            end
        end
    end   
end

function testbackup_w!(f::Fluid)
    @sync for pid in workers()
        @spawnat pid begin
            @inbounds @simd for c in localpart(f.cells)
                c.wb = copy(c.w)
            end
        end
    end   
end

function copyfluid!(f::Fluid)
    f1 = Fluid(f.dim)
    f1.point1 = copy(f.point1)
    f1.point2 = copy(f.point2)
    f1.nmesh = copy(f.nmesh)
    f1.d = copy(f.d)
    f1.ng = copy(f.ng)
    f1.dist = copy(f.dist)
    f1.ndiv = copy(f.ndiv)
    f1.cells = copy!(f.cells, f.ndiv, f.dist)
    f1.boundaries = copy(f.boundaries)
    f1.para = f.para
    return f1
end

function Base.copy!(c::Cell) 
    cell = Cell(length(c.u))
    # 必须复制坐标
    cell.rho = c.rho
    cell.u = copy(c.u)
    cell.e = c.e
    cell.p = c.p
    cell.w = copy(c.w)
    cell.wb = copy(c.wb)   
    cell.rhs = copy(c.rhs)
    cell.i = copy(c.i)
    cell.x = copy(c.x)
    cell.mark = c.mark
    cell.target_id = copy(c.target_id)
    return cell
end


function Base.copy!(a::Array{Vector{Union{Nothing, Float64}}})
    b = Array{Vector{Union{Nothing, Float64}}}(undef, size(a))
    for i in eachindex(a)
        b[i] = copy(a[i])
    end
    return b
end

function Base.copy!(a::Vector{Union{Nothing, Float64}})
    b = Vector{Union{Nothing, Float64}}(undef, size(a))
    for i in eachindex(a)
        b[i] = a[i]
    end
    return b
end

function Base.copy!(a::Array{Cell}) 
    b = Array{Cell}(undef, size(a))
    for k in eachindex(b)
        b[k] = copy!(a[k])
    end
    return b
end

"""
这个函数如果写的方式不对，就会相当耗时。
"""
function Base.copy!(a::DArray{Cell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
    dim = length(a[1].u)
    tmp = Array{Cell}(undef, size(a))
    for i in eachindex(tmp)
        tmp[i] = Cell(dim)
    end
    b = distribute(tmp, procs = workers(), dist = dist)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(a)
                copy_complete_cell_to!(c, b[Tuple(c.i)...])
            end
        end
    end
    return b
end

"""
这个函数如果写的方式不对，就会相当耗时。
"""
function testcopy!(a::DArray{Cell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
    dim = a[1].dim
    tmp = Array{Cell}(undef, size(a))
    for i in eachindex(tmp)
        tmp[i] = Cell(dim)
    end
    b = distribute(tmp, procs = workers(), dist = dist)
    if dim == 1
        @sync for pid in workers()
            @spawnat pid begin
                for c in localpart(a)
                    copy_complete_cell_to!(c, b[Tuple(c.i)...])
                end 
            end
        end
    elseif dim == 2
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(a)
                for i in inds[1]
                    for j in inds[2]
                        copy_complete_cell_to!(a[i, j], b[i, j])
                    end
                end 
            end
        end
    elseif dim == 3
        @sync for pid in workers()
            @spawnat pid begin
                inds = localindices(a)
                for i in inds[1]
                    for j in inds[2]
                        for k in inds[3]
                            copy_complete_cell_to!(a[i, j, k], b[i, j, k])
                        end
                    end
                end 
            end
        end
    else
        error("undef dim")
    end
    return b
end

# function Base.copy!(a::DArray{Cell, N, A} where N where A, ndiv::Int, dist::Array{Int}) 
#     tmp = Array{Cell}(undef, size(a))
#     @time for i in eachindex(a)
#         c = a[i]
#         tmp[i] = copy!(c)
#     end
    
#     b = distribute(tmp, procs = workers(), dist = dist)
#     return b
# end

function Base.copy!(a::SubArray{Cell,N,P,I,L} where L where I where P<:DArray where N)
    b = Array{Cell}(undef, size(a))
    for k in eachindex(b)
        b[k] = copy!(a[k])
    end
    return b
end

function copy_to!(c::Cell, target::Cell)
    target.rho = c.rho
    target.u = copy(c.u)
    target.e = c.e
    target.p = c.p
    target.w = copy(c.w)
    target.wb = copy(c.wb)
end

function copy_to!(a::Array{Cell}, target::Array{Cell}) 
    @assert size(a) == size(target)
    for k in eachindex(target)
        copy_to!(a[k], target[k])
    end
end
function copy_to!(a::SubArray, target::SubArray) 
    @assert size(a) == size(target)
    for k in eachindex(target)
        copy_to!(a[k], target[k])
    end
end
function copy_to!(a::SubArray, target::Array{Cell}) 
    @assert size(a) == size(target)
    for k in CartesianIndices(target)
        # println("-- copy_to: 1 --")
        # println(k)
        # println(typeof(a))
        # println(typeof(target))
        # println(typeof(a[k]))
        # println(typeof(target[k]))
        copy_to!(a[k], target[k])
    end
end

function copy_complete_cell_to!(c::Cell, target::Cell)
    target.rho = c.rho
    target.u = copy(c.u)
    target.e = c.e
    target.p = c.p
    target.w = copy(c.w)
    target.wb = copy(c.wb)
    target.rhs = copy(c.rhs)
    target.i = copy(c.i)
    target.x = copy(c.x)
end