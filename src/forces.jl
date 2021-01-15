using ForwardDiff
using Setfield
using StaticArrays

include("energies.jl")



""" Force between two atoms"""
@inline function pair_force(r1::SVector{3, Float64}, r2::SVector{3, Float64}, box_size)
    """
    Calculates the Lennard Jones Force between two atoms

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
    f1 : SVector{3}
        SVector with x, y, z forces on atom 1
    f2 : SVector{3}
        SVector with x, y, z forces on atom 2

    """
    diff = SVector{3}(0.0, 0.0, 0.0)

    # apply mirror image separation
    @inbounds for i=1:3
            diff = @set diff[i] = vector1D(r1[i], r2[i], box_size[i])
        end

    rij_sq = diff[1]*diff[1] + diff[2]*diff[2] + diff[3]*diff[3]

    if  rij_sq > box_size[1] / 2
        return SVector{3}(0.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0)
    else

        sr2 = 1 / rij_sq #σ^2 / rij_sq
        sr6 = sr2 ^ 3
        sr12 = sr6 ^ 2
        f = (24 / rij_sq) * (2 * sr12 - sr6) # ((24 * ϵ) / rij_sq) * (2 * sr12 - sr6)

        dx = -f * diff[1]
        dy = -f * diff[2]
        dz = -f * diff[3]

        f1 = SVector{3}(dx, dy, dz)
        f2 = SVector{3}(-dx, -dy, -dz)
    end

    # return forces between atom 1 and atom 2
    return f1, f2
end


function numerical_total_force(r::Vector{SVector{3, Float64}}, box_size::SVector{3,Float64})
    """
    Calculates the Lennard Jones forces on all atoms using numerical derivatives

    Parameters
    ----------
    r : SVector{3}
        coordinates of all atoms
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    forces : Vector{SVector{3}}
        Vector of forces, where each element is the x, y, z forces on that atom
    """
    n = length(r)
    dx = 1e-8
    forces = [SVector{3}(0.0, 0.0, 0.0) for i=1:n ]
    for i = 1:(n-1)
        for j = (i+1):n
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[1] += dx
            rxm = @set rxm[1] -= dx
            dE_dx = -( pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size) ) / (2 * dx)
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[2] += dx
            rxm = @set rxm[2] -= dx
            dE_dy = -( pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size) ) / (2 * dx)
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[3] += dx
            rxm = @set rxm[3] -= dx
            dE_dz = -( pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size) ) / (2 * dx)
            forces[i] = forces[i] + SVector{3}(dE_dx, dE_dy, dE_dz)
            forces[j] = forces[j] - SVector{3}(dE_dx, dE_dy, dE_dz)
        end
    end
    return forces
end

# analytical_force = x -> ForwardDiff.gradient(pair_energy, x, y, box_size)
""" Calculates the force between two atoms using ForwardDiff"""
grad(x,y,b) = -ForwardDiff.gradient(x -> pair_energy(x, y, b), x)

function analytical_total_force(r::Vector{SVector{3, Float64}}, box_size::SVector{3,Float64})
    """
    Calculates the Lennard Jones forces on all atoms using AutoDifferentiation

    Parameters
    ----------
    r : SVector{3}
        coordinates of all atoms
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    forces : Vector{SVector{3}}
        Vector of forces, where each element is the x, y, z forces on that atom
    """
    n = length(r)
    forces = [SVector{3}(0.0, 0.0, 0.0) for i=1:n ]

    for i = 1:(n-1)
        for j = (i+1):n
            dE_dr = grad(r[i], r[j], box_size)  #
            forces[i] = forces[i] + dE_dr
            forces[j] = forces[j] - dE_dr
        end
    end
    return forces
end

function total_force(r::Vector{SVector{3, Float64}}, box_size::SVector{3,Float64})
    """
    Calculates the Lennard Jones forces on all atoms using by-hand
    analytical derivatives

    Parameters
    ----------
    r : SVector{3}
        coordinates of all atoms
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    forces : Vector{SVector{3}}
        Vector of forces, where each element is the x, y, z forces on that atom
    """
    n = length(r)
    forces = [SVector{3}(0.0, 0.0, 0.0) for i=1:n ]
    for i = 1:(n-1)
        for j = (i+1):n
            fi, fj = pair_force(r[i], r[j], box_size)
            forces[i] = forces[i] + fi
            forces[j] = forces[j] + fj
        end
    end
    return forces
end
