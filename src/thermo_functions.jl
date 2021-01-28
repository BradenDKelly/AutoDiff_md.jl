#using StaticArrays
#using Statistics
export
    virial,
    press_full,
    kinetic_energy

include("forces.jl")

"Vector between two coordinate values, accounting for mirror image seperation"
@inline function vector1D(c1, c2, box_size)
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end
function virial(r::Vector, eps::Vector, sig::Vector, cutoff::Real, box_size)
    """ Calculates the virial contribution to pressure

    Parameters
    ----------
    r : Vector{SVector{3}}
        Vector of atom coordinates
    eps : Vector
        epsilon parameters of all atoms
    sigma : Vector
        sigma parameters of all atoms
    cutoff : Float
        interaction cutoff distance
    box : SVector{3}
        Vector with x, y, z box lengths

    Returns
    ----------
    tot_vir : Float64
        total virial contribution to the pressure
    """
    n = length(r)
    tot_vir = 0.0
    for j=1:n-1
        for k=j+1:n
            fij = grad(r[k],r[j], eps[j], eps[k], sig[j], sig[k], cutoff, box_size)
            dx = vector1D(r[j][1], r[k][1], box_size[1])
            dy = vector1D(r[j][2], r[k][2], box_size[2])
            dz = vector1D(r[j][3], r[k][3], box_size[3])
            #println(i, j, fij, dx, dy, dz)
            vir = dot(fij, SVector([dx, dy, dz]...) )
            tot_vir += vir
        end
    end
    return tot_vir
end
press_full(vir::Float64, n::Int64, vol::Float64, tmp::Float64) = n * 0.0083144621 / vol * tmp + vir / vol / 3
kinetic_energy(v::Vector, m:: Vector) = 0.5 * dot(m .* v,v)
