#using StaticArrays
#using Statistics
export
    virial,
    press_full,
    kinetic_energy,
    volume

include("forces.jl")
include("structs.jl")


@inline function vector1D(c1, c2, box_size)
    """Vector between two coordinate values, accounting for mirror image seperation"""
    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) : (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) : (c2 - c1 + box_size)
    end
end

function virial(simulation_arrays::SimulationArrays, cutoff, box_size)
    return lj_virial(simulation_arrays, cutoff, box_size)
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

function lj_virial(simulation_arrays::SimulationArrays, cutoff::Real, box_size)
    """ Calculates the virial contribution to pressure

    Parameters
    ----------
    simulation_arrays : SimulationArrays
        atomic and molecular arrays with r, v, m, f, etc.
    cutoff : Float
        interaction cutoff distance
    box : SVector{3}
        Vector with x, y, z box lengths

    Returns
    ----------
    tot_vir : Float64
        total molecular virial contribution to the pressure
    """
    n = length(simulation_arrays.molecule_arrays.COM[:])
    tot_vir = 0.0
    for i=1:n-1
        atomIs = simulation_arrays.molecule_arrays[i].firstAtom
        atomIe = simulation_arrays.molecule_arrays[i].lastAtom
        ri = Center_of_Mass(
            simulation_arrays.atom_arrays.r[atomIs:atomIe],
            simulation_arrays.atom_arrays.mass[atomIs:atomIe]
            )
        for j=i+1:n
            atomJs = simulation_arrays.molecule_arrays[j].firstAtom
            atomJe = simulation_arrays.molecule_arrays[j].lastAtom
            rj = Center_of_Mass(
                simulation_arrays.atom_arrays.r[atomJs:atomJe],
                simulation_arrays.atom_arrays.mass[atomJs:atomJe]
                )
            dx = vector1D(ri[1], rj[1], box_size[1])
            dy = vector1D(ri[2], rj[2], box_size[2])
            dz = vector1D(ri[3], rj[3], box_size[3])
            rij = SVector(dx, dy, dz)
            # cycle over atoms in first molecule
            for a = atomIs:atomIe
                ta = simulation_arrays.atom_arrays.atype[a]
                # cycle over atoms in second molecule
                for b = atomJs:atomJe
                    tb = simulation_arrays.atom_arrays.atype[b]
                    fab = lj_grad(
                        simulation_arrays.atom_arrays.r[a],
                        simulation_arrays.atom_arrays.r[b],
                        simulation_arrays.vdwTable.ϵᵢⱼ[ta, tb],
                        simulation_arrays.vdwTable.σᵢⱼ[ta, tb],
                        cutoff,
                        box_size
                        )
                    # project atom-atom force onto molecule-molecule vector
                    vir = dot(fab, rij)
                    tot_vir += vir
                end # for b
            end # for a
        end # for j
    end # for i
    return tot_vir
end

volume(box_size::SVector) = box_size[1] * box_size[2] * box_size[3]
press_full(vir::Float64, n::Int64, vol::Float64, tmp::Float64) = n * 0.0083144621 / vol * tmp + vir / vol / 3
kinetic_energy(v::Vector, m:: Vector) = 0.5 * dot(m .* v,v)
