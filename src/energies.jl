


""" Potential energy between two atoms"""
@inline function pair_energy(r1::SVector{3}, r2::SVector{3}, box_size::SVector{3})
    """
    Calculates the Lennard Jones energy between two atoms

    Parameters
    ----------
    r1 : SVector{3}
        coordinate of atom 1
    r2 : SVector{3}
        coordinate of atom 2
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    e : float
        potential energy between atoms 1 and atom 2
    """
    # diff = SVector{3}(0.0,0.0,0.0)
    # apply mirror image separation
    #@inbounds for i=1:
    #        diff = @set diff[i] = vector1D(r1[i], r2[i], box_size[i])
    #end
    dx = vector1D(r1[1], r2[1], box_size[1])
    dy = vector1D(r1[2], r2[2], box_size[2])
    dz = vector1D(r1[3], r2[3], box_size[3])

    rij_sq = dx * dx + dy * dy + dz * dz
    #rij_sq = diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3]

    if  rij_sq > 2.5^2 #rij_sq > box_size[1] / 2
        return 0.0 * rij_sq
    end

    sr2 = 1 / rij_sq
    sr6 = sr2^3
    sr12 = sr6^2
    # to match allen & Tildesly I add the cut potential manually
    # do not actually need 4 * (-0.004079222784000001)  otherwise
    e = 4 * (sr12 - sr6) - 4 * (-0.004079222784000001)
    return e

end

using StaticArrays

function total_energy(r::Vector{SVector{3, Float64}}, box_size::SVector{3,Float64})
    """
    Calculates the Lennard Jones energy of the system

    Parameters
    ----------
    r : SVector{3}
        coordinates of all atoms
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    energy : float
        LJ potential energy of the system
    """
    n = length(r)
    energetics = [0.0 for i=1:n ]
    energy = 0.0
    for i = 1:(n-1)
        for j = (i+1):n
            ans = pair_energy(r[i], r[j], box_size)
            energy += ans
            energetics[i] += ans
        end
    end
    return energy, energetics
end
