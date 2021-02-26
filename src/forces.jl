export
    bond_force,
    bond_grad,
    total_force_bond,
    grad,
    lj_grad,
    total_lj_force,
    pair_force,
    numerical_total_force,
    analytical_total_force,
    total_force,
    vector

include("energies.jl")

#""" Calculates the force between two atoms using ForwardDiff"""
grad(x, y, e1, e2, s1, s2, c, b) =
    -ForwardDiff.gradient(x -> pair_energy(x, y, e1, e2, s1, s2, c, b), x)

#""" Calculates the force between two atoms using ForwardDiff"""
lj_grad(x, y, e, s, c, b) =
    -ForwardDiff.gradient(x -> lj_atom_pair_energy(x, y, e, s, c, b), x)

#""" Calculates the force between two atoms using ForwardDiff"""
lj_grad(x, y, e, s, c, b, shift) =
    -ForwardDiff.gradient(x -> lj_atom_pair_energy(x, y, e, s, c, b, shift), x)

bond_grad(x, y, k, l, b) =
    -ForwardDiff.gradient(x -> spring_energy(x, y, k, l, b), x)

#"3D vector between two `Coordinates`, accounting for mirror image sep."
vector(coords_one::SVector, coords_two::SVector, box_size::SVector) = [
    vector1D(coords_one[1], coords_two[1], box_size[1]),
    vector1D(coords_one[2], coords_two[2], box_size[2]),
    vector1D(coords_one[3], coords_two[3], box_size[3]),
]

#"Force on each atom in a covalent bond."
function bond_force(
    coords_one::SVector,
    coords_two::SVector,
    springConstant::F,
    bondLength::F,
    box_size::SVector,
) where F

    ab = vector(coords_one, coords_two, box_size)
    c = springConstant * (norm(ab) - bondLength)
    f1 = c * normalize(ab)

    return f1, -f1
end

#"""Calculates the forces on all atoms due to covalent bonds"""
function total_force_bond(
    sa::SimulationArrays,
    box_size::SVector,
)
    n = length(sa.atom_arrays.r[:])
    fa = zero(Vector{Float64}(undef, 3))
    fb = zero(Vector{Float64}(undef, 3))
    bond_forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    for i = 1:length(sa.intraFF.bonds)
        ai = sa.intraFF.bonds[i].ai
        aj = sa.intraFF.bonds[i].aj

        #fa, fb = bond_force(
        fa = bond_grad(
            sa.atom_arrays[ai].r,
            sa.atom_arrays[aj].r,
            sa.intraFF.bonds[i].kparam,
            sa.intraFF.bonds[i].bondLength,
            box_size,
        )
        bond_forces[ai] += fa
        #bond_forces[aj] += fb
        bond_forces[aj] += -fa
        # atomArray[ai].f.x += fa[1]
        # atomArray[aj].f.x += fb[1]
        # atomArray[ai].f.y += fa[2]
        # atomArray[aj].f.y += fb[2]
        # atomArray[ai].f.z += fa[3]
        # atomArray[aj].f.z += fb[3]

    end
    return bond_forces
end
#=
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
=#
@inline function total_lj_force2(
    simulation_arrays::SimulationArrays,
    cutoff::T,
    box_size::SVector,
) where {T}
    n::Int64 = length(simulation_arrays.atom_arrays.r[:])
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    ti::Int64 = 0
    tj::Int64 = 0

    for i = 1:(n-1)
        ti = simulation_arrays.atom_arrays.atype[i]
        for j = (i+1):n
            tj = simulation_arrays.atom_arrays.atype[j]
            dE_dr = lj_grad(
                simulation_arrays.atom_arrays.r[i],
                simulation_arrays.atom_arrays.r[j],
                simulation_arrays.vdwTable.ϵᵢⱼ[ti, tj],
                simulation_arrays.vdwTable.σᵢⱼ[ti, tj],
                cutoff,
                box_size,
            )
            forces[i] = forces[i] + dE_dr
            forces[j] = forces[j] - dE_dr
        end
    end
    return forces
end
#=
"""
Calculates the Lennard Jones forces on all atoms using AutoDifferentiation
and neighborlist

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
=#
@inline function total_lj_force(
    sa::SimulationArrays,
    cutoff::T,
    box_size::SVector,
    point::Array,
    list::Array,
) where {T}
    n = length(sa.atom_arrays.r[:])
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    for i = 1:(n-1)
        ti = sa.atom_arrays.atype[i]
        for j = point[i]:(point[i+1]-1)
            k = list[j]
            tj = sa.atom_arrays.atype[k]
            dE_dr = lj_grad(
                sa.atom_arrays.r[i],
                sa.atom_arrays.r[k],
                sa.vdwTable.ϵᵢⱼ[ti, tj],
                sa.vdwTable.σᵢⱼ[ti, tj],
                cutoff,
                box_size,
            )
            forces[i] = forces[i] + dE_dr
            forces[k] = forces[k] - dE_dr
        end
    end
    return forces
end
#=
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
=#
@inline function total_lj_force(
    r::Vector,
    atype::Vector,
    vdwTable::Tables,
    cutoff::T,
    box_size::SVector,
) where {T}

    n::Int64 = length(r)
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    ti::Int64 = 0
    tj::Int64 = 0

    for i = 1:(n-1)
        ti = atype[i]
        for j = (i+1):n
            tj = atype[j]
            dE_dr = lj_grad(
                r[i],
                r[j],
                vdwTable.ϵᵢⱼ[ti, tj],
                vdwTable.σᵢⱼ[ti, tj],
                cutoff,
                box_size,
            )
            forces[i] = forces[i] + dE_dr
            forces[j] = forces[j] - dE_dr
        end
    end
    return forces
end

#"""total energy of the system"""
function analytical_total_force(
    simulation_arrays::SimulationArrays,
    cutoff,
    box_size,
    point,
    list,
)
    bond_forces = total_force_bond(simulation_arrays, box_size)
    #println("bonds: ",bond_forces[1])
    lj_forces = total_lj_force(simulation_arrays, cutoff, box_size, point, list)
    #println("lj: ",lj_forces[1])
    return bond_forces + lj_forces
end

#"""total energy of the system"""
function analytical_total_force(
    r::Vector,
    atype::Vector,
    vdwTable::Tables,
    cutoff,
    box_size,
)

    return total_lj_force(r, atype, vdwTable, cutoff, box_size)
end
#=
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
=#
@inline function pair_force(
    r1::SVector{3,Float64},
    r2::SVector{3,Float64},
    box_size,
)
    diff = SVector{3}(0.0, 0.0, 0.0)

    # apply mirror image separation
    @inbounds for i = 1:3
        diff = @set diff[i] = vector1D(r1[i], r2[i], box_size[i])
    end

    rij_sq = diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3]

    if rij_sq > box_size[1] / 2
        return SVector{3}(0.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0)
    else

        sr2 = 1 / rij_sq #σ^2 / rij_sq
        sr6 = sr2^3
        sr12 = sr6^2
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
#=
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
=#
function numerical_total_force(
    r::Vector{SVector{3,Float64}},
    box_size::SVector{3,Float64},
)
    n = length(r)
    dx = 1e-8
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[1] += dx
            rxm = @set rxm[1] -= dx
            dE_dx =
                -(
                    pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size)
                ) / (2 * dx)
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[2] += dx
            rxm = @set rxm[2] -= dx
            dE_dy =
                -(
                    pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size)
                ) / (2 * dx)
            rxp, rxm = deepcopy(r[i]), deepcopy(r[i])
            rxp = @set rxp[3] += dx
            rxm = @set rxm[3] -= dx
            dE_dz =
                -(
                    pair_energy(rxp, r[j], box_size) -
                    pair_energy(rxm, r[j], box_size)
                ) / (2 * dx)
            forces[i] = forces[i] + SVector{3}(dE_dx, dE_dy, dE_dz)
            forces[j] = forces[j] - SVector{3}(dE_dx, dE_dy, dE_dz)
        end
    end
    return forces
end
#=
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
=#
function analytical_total_force(
    r::Vector{SVector{3,Float64}},
    eps::Vector,
    sig::Vector,
    cutoff::Real,
    box_size::SVector{3,Float64},
)

    n = length(r)
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]

    for i = 1:(n-1)
        for j = (i+1):n
            dE_dr = grad(
                r[i],
                r[j],
                eps[i],
                eps[j],
                sig[i],
                sig[j],
                cutoff,
                box_size,
            )  #
            forces[i] = forces[i] + dE_dr
            forces[j] = forces[j] - dE_dr
        end
    end
    return forces
end

#=
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
=#
function total_force(
    r::Vector{SVector{3,Float64}},
    box_size::SVector{3,Float64},
)

    n = length(r)
    forces = [SVector{3}(0.0, 0.0, 0.0) for i = 1:n]
    for i = 1:(n-1)
        for j = (i+1):n
            fi, fj = pair_force(r[i], r[j], box_size)
            forces[i] = forces[i] + fi
            forces[j] = forces[j] + fj
        end
    end
    return forces
end
