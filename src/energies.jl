export pair_energy, total_energy

include("structs.jl")

@inline function lj_atom_pair_energy(
    r1::SVector{3},
    r2::SVector{3},
    ϵ::Real,
    σ::Real,
    cutoff::Real,
    box_size::SVector{3},
    shifted_potential::Real = 0.0,
)
    """
    Calculates the Lennard Jones energy between two atoms

    Parameters
    ----------
    r1 : SVector{3}
        coordinate of atom 1
    r2 : SVector{3}
        coordinate of atom 2
    ϵ : Scalar
        post-mixing rule epsilon parameter of atom pair
    σ : Vector
        post-mixing rule sigma parameter of atom pair
    cutoff : Float
        interaction cutoff distance
    box_size : SVector{3}
        vector with box length in the x, y, z direction
    shifted_potential : Float
        cut and shift potential

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

    # TODO precalculate squared cutoff
    if rij_sq > cutoff^2 #rij_sq > box_size[1] / 2
        return 0.0 * rij_sq
    end

    # TODO precalculate squared sigma
    sr2 = σ^2 / rij_sq
    sr6 = sr2^3
    sr12 = sr6^2
    # to match allen & Tildesly I add the cut potential manually
    # do not actually need 4 * (-0.004079222784000001)  otherwise
    # current shift is for r=1.2 nm
    e = 4 * ϵ * (sr12 - sr6) - shifted_potential #(-0.002078435714914992) #- 4 * (-0.004079222784000001)
    return e
end

@inline function lj_molec_pair_energy()
    println("In progress")
end

@inline function lj_molec_vs_system_energy()
    println("In progress")

end

@inline function total_lj_energy(
    simulation_arrays::SimulationArrays,
    cutoff::T,
    box_size::SVector,
) where {T}
    """
    Calculates the Lennard Jones forces on all atoms using AutoDifferentiation

    Parameters
    ----------
    simulation_array : SimulationArray
        atom_arrays::StructArray
            molNum::I
            molType::I
            atype::I
            mass::F
            r::SVector{3,F}
            v::SVector{3,F}
            f::SVector{3,F}
            qq::Float64
    cutoff : Float
        interaction cutoff
    box_size : SVector{3}
        vector with box length in the x, y, z direction

    Returns
    ---------
    forces : Vector{SVector{3}}
        Vector of forces, where each element is the x, y, z forces on that atom
    """
    n = length(simulation_arrays.atom_arrays.r[:])
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    energetics = [0.0 for i = 1:n]
    energy = 0.0
    for i = 1:(n-1)
        ti = simulation_arrays.atom_arrays.atype[i]
        for j = (i+1):n
            tj = simulation_arrays.atom_arrays.atype[j]
            ener = lj_atom_pair_energy(
                simulation_arrays.atom_arrays.r[i],
                simulation_arrays.atom_arrays.r[j],
                simulation_arrays.vdwTable.ϵᵢⱼ[ti, tj],
                simulation_arrays.vdwTable.σᵢⱼ[ti, tj],
                cutoff,
                box_size,
            )
            energy += ener
            energetics[i] += ener
        end
    end
    return energy, energetics
end

function total_energy(simulation_arrays::SimulationArrays, cutoff, box_size)
    """total energy of the system"""
    return total_lj_energy(simulation_arrays, cutoff, box_size)
end

@inline function pair_energy(
    r1::SVector{3},
    r2::SVector{3},
    eps1::Real,
    eps2::Real,
    sig1::Real,
    sig2::Real,
    cutoff::Real,
    box_size::SVector{3},
)
    """
    NOTE: OBSOLETE FUNCTION:

    Calculates the Lennard Jones energy between two atoms

    Parameters
    ----------
    r1 : SVector{3}
        coordinate of atom 1
    r2 : SVector{3}
        coordinate of atom 2
    eps : Vector
        epsilon parameters of all atoms
    sigma : Vector
        sigma parameters of all atoms
    cutoff : Float
        interaction cutoff distance
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

    if rij_sq > cutoff^2 #rij_sq > box_size[1] / 2
        return 0.0 * rij_sq
    end
    sigma = (sig1 + sig2) / 2
    epsilon = sqrt(eps1 * eps2)

    sr2 = sigma^2 / rij_sq
    sr6 = sr2^3
    sr12 = sr6^2
    # to match allen & Tildesly I add the cut potential manually
    # do not actually need 4 * (-0.004079222784000001)  otherwise
    # current shift is for r=1.2 nm
    e = 4 * epsilon * (sr12 - sr6) - (-0.002078435714914992) #- 4 * (-0.004079222784000001)
    return e

end

#using StaticArrays

function total_energy(
    r::Vector{SVector{3,Float64}},
    eps::Vector,
    sig::Vector,
    cutoff::Real,
    box_size::SVector{3,Float64},
)
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
    energetics = [0.0 for i = 1:n]
    energy = 0.0
    @inbounds for i = 1:(n-1)
        @inbounds @simd for j = (i+1):n
            ans = pair_energy(
                r[i],
                r[j],
                eps[i],
                eps[j],
                sig[i],
                sig[j],
                cutoff,
                box_size,
            )
            energy += ans
            energetics[i] += ans
        end
    end
    return energy, energetics
end
