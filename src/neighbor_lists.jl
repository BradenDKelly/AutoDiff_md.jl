export VerletList, make_neighbor_list!

mutable struct VerletList{T} <: NeighborList
    buffer::T
    point::Array
    list::Array
end

# TODO VerletList{Any} should be fixed to include a type T
function make_neighbor_list!(
    r::Vector,
    nonbonded_matrix::Array{Bool,2},
    box_size::SVector{3,T},
    cutoff::T,
    array_holder::VerletList{T},
) where {T}
    """
    Generates a verlet neighbor list.

    Fortran Code https://github.com/Allen-Tildesley/examples/blob/master/verlet_list_module.f90
    # They also have python code too
    # In: r,r_list_box, BoxSize, nonbonded_matrix
    # Out: list,point
    """

    k = 0
    n = length(r)
    point = zeros(MVector{n,Int64})
    diff = MVector{3}(0.0, 0.0, 0.0)
    r_list_box_sq = (cutoff + array_holder.buffer)^2
    mlist = Int64[]

    for i = 1:(n-1) # ! Begin outer loop over atoms

        point[i] = k + 1

        for j = (i+1):n #! Begin inner loop over partner atoms

            @inbounds for l = 1:3
                diff[l] = vector1D(r[i][l], r[j][l], box_size[l])
            end
            rij_sq = diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3]

            if (rij_sq < r_list_box_sq && nonbonded_matrix[i, j])
                k = k + 1
                push!(mlist, j)

            end

        end #! End inner loop over partner atoms

    end #! End outer loop over atoms
    list = Vector{Int64}(undef, length(mlist))
    #list = zeros(MVector{length(mlist), Int64})
    #list = MVector{length(mlist)}(mlist...)
    for i = 1:length(mlist)
        list[i] = mlist[i]
    end
    point[n] = k + 1
    array_holder.point = point
    array_holder.list = list
    # = VerletList(array_holder.buffer, point, list)
    #return VerletList(array_holder.buffer, point, list)
end
