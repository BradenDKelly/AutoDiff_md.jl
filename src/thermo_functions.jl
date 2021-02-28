export virial,
    press_full,
    kinetic_energy,
    volume,
    lj_virial,
    press_full_tensor,
    lj_virial_tensor,
    kinetic_energy_tensor,
    pressure_from_tensors

include("forces.jl")
include("structs.jl")

# TODO hand code LJ force calc... it is 2-3x faster than autodiff :(
#"""Vector between two coordinate values, accounting for mirror image seperation"""
@inline function vector1D(c1, c2, box_size)

    if c1 < c2
        return (c2 - c1) < (c1 - c2 + box_size) ? (c2 - c1) :
               (c2 - c1 - box_size)
    else
        return (c1 - c2) < (c2 - c1 + box_size) ? (c2 - c1) :
               (c2 - c1 + box_size)
    end
end

function virial(simulation_arrays::SimulationArrays, cutoff, box_size)
    return lj_virial(simulation_arrays, cutoff, box_size)
end
#=
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

NOTES: Allen & Tildesley 2nd Ed. pp. 62 for a good reference
Molecular systems can be calculated using atom-atom forces - BUT - this
requires using TOTAL force i.e., include force due to bonds/angles/torsions

optionally - can use molecular version and project atom-atom forces onto
molecule-molecule rij

If using constraints, must include constraint contribution to virial
"""
=#
function virial(r::Vector, eps::Vector, sig::Vector, cutoff::Real, box_size)

    n = length(r)
    tot_vir = 0.0
    @inbounds for j = 1:n-1
        @inbounds for k = j+1:n
            fij = grad(
                r[k],
                r[j],
                eps[j],
                eps[k],
                sig[j],
                sig[k],
                cutoff,
                box_size,
            )
            dx = vector1D(r[j][1], r[k][1], box_size[1])
            dy = vector1D(r[j][2], r[k][2], box_size[2])
            dz = vector1D(r[j][3], r[k][3], box_size[3])
            #println(i, j, fij, dx, dy, dz)
            vir = dot(fij, SVector([dx, dy, dz]...))
            tot_vir += vir
        end
    end
    return tot_vir
end

#=
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

    NOTES: Allen & Tildesley 2nd Ed. pp. 62 for a good reference
    Molecular systems can be calculated using atom-atom forces - BUT - this
    requires using TOTAL force i.e., include force due to bonds/angles/torsions

    optionally - can use molecular version and project atom-atom forces onto
    molecule-molecule rij

    If using constraints, must include constraint contribution to virial
"""
=#
function lj_virial(simulation_arrays::SimulationArrays, cutoff::Real, box_size)

    n = length(simulation_arrays.molecule_arrays.COM[:])
    tot_vir = 0.0
    @inbounds for i = 1:n-1
        atomIs = simulation_arrays.molecule_arrays[i].firstAtom
        atomIe = simulation_arrays.molecule_arrays[i].lastAtom
        ri = Center_of_Mass(
            simulation_arrays.atom_arrays.r[atomIs:atomIe],
            simulation_arrays.atom_arrays.mass[atomIs:atomIe],
        )
        @inbounds for j = i+1:n
            atomJs = simulation_arrays.molecule_arrays[j].firstAtom
            atomJe = simulation_arrays.molecule_arrays[j].lastAtom
            rj = Center_of_Mass(
                simulation_arrays.atom_arrays.r[atomJs:atomJe],
                simulation_arrays.atom_arrays.mass[atomJs:atomJe],
            )
            dx = vector1D(ri[1], rj[1], box_size[1])
            dy = vector1D(ri[2], rj[2], box_size[2])
            dz = vector1D(ri[3], rj[3], box_size[3])
            rij = SVector(dx, dy, dz)
            #if norm(rij) < cutoff
            # cycle over atoms in first molecule
            mol_vir = 0.0
            @inbounds for a = atomIs:atomIe
                ta = simulation_arrays.atom_arrays.atype[a]
                # cycle over atoms in second molecule
                @inbounds for b = atomJs:atomJe
                    tb = simulation_arrays.atom_arrays.atype[b]
                    fab = lj_grad(
                        simulation_arrays.atom_arrays.r[a],
                        simulation_arrays.atom_arrays.r[b],
                        simulation_arrays.vdwTable.ϵᵢⱼ[ta, tb],
                        simulation_arrays.vdwTable.σᵢⱼ[ta, tb],
                        cutoff,
                        box_size,
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

function lj_virial_tensor(
    simulation_arrays::SimulationArrays,
    cutoff::Real,
    box_size,
)

    n = length(simulation_arrays.molecule_arrays.COM[:])

    mol_vir = 0.0
    virial_tensor = zeros(MMatrix{3,3})
    @inbounds for i = 1:n-1
        atomIs = simulation_arrays.molecule_arrays[i].firstAtom
        atomIe = simulation_arrays.molecule_arrays[i].lastAtom
        ri = Center_of_Mass(
            simulation_arrays.atom_arrays.r[atomIs:atomIe],
            simulation_arrays.atom_arrays.mass[atomIs:atomIe],
        )
        @inbounds for j = i+1:n
            atomJs = simulation_arrays.molecule_arrays[j].firstAtom
            atomJe = simulation_arrays.molecule_arrays[j].lastAtom
            rj = Center_of_Mass(
                simulation_arrays.atom_arrays.r[atomJs:atomJe],
                simulation_arrays.atom_arrays.mass[atomJs:atomJe],
            )
            dx = vector1D(ri[1], rj[1], box_size[1])
            dy = vector1D(ri[2], rj[2], box_size[2])
            dz = vector1D(ri[3], rj[3], box_size[3])
            rij = SVector(dx, dy, dz)
            #if norm(rij) < cutoff
            # cycle over atoms in first molecule
            mol_vir = 0.0
            @inbounds for a = atomIs:atomIe
                ta = simulation_arrays.atom_arrays.atype[a]
                ra = simulation_arrays.atom_arrays.r[a]
                # cycle over atoms in second molecule
                @inbounds for b = atomJs:atomJe
                    tb = simulation_arrays.atom_arrays.atype[b]
                    rb = simulation_arrays.atom_arrays.r[b]
                    fab = lj_grad(
                        simulation_arrays.atom_arrays.r[a],
                        simulation_arrays.atom_arrays.r[b],
                        simulation_arrays.vdwTable.ϵᵢⱼ[ta, tb],
                        simulation_arrays.vdwTable.σᵢⱼ[ta, tb],
                        cutoff,
                        box_size,
                    )
                    # project atom-atom force onto molecule-molecule vector
                    vir = dot(fab, rij)
                    mol_vir += vir

                    dx = vector1D(ra[1], rb[1], box_size[1])
                    dy = vector1D(ra[2], rb[2], box_size[2])
                    dz = vector1D(ra[3], rb[3], box_size[3])
                    rab = SVector(dx, dy, dz)

                    virial_tensor += rab * fab' # outer product, 3x3 matrix
                end # for b
            end # for a
        end # for j
    end # for i
    #println("virial tensor: $virial_tensor")
    return virial_tensor
end

function total_virial_bond_tensor(sa::SimulationArrays, box_size::SVector)
    n = length(sa.atom_arrays.r[:])
    fab = zero(Vector{Float64}(undef, 3))

    virial_tensor_bond = zeros(MMatrix{3,3})
    @inbounds for i = 1:length(sa.intraFF.bonds)
        ai = sa.intraFF.bonds[i].ai
        aj = sa.intraFF.bonds[i].aj

        fab = bond_grad(
            sa.atom_arrays[ai].r,
            sa.atom_arrays[aj].r,
            sa.intraFF.bonds[i].kparam,
            sa.intraFF.bonds[i].bondLength,
            box_size,
        )
        rab = vector(sa.atom_arrays[ai].r, sa.atom_arrays[aj].r, box_size)
        virial_tensor_bond += rab * fab'   # outer product 3x3 matrix
    end
    return virial_tensor_bond
end


function kinetic_energy_tensor(sa::SimulationArrays)
    kinetic_energy = zeros(MMatrix{3,3})
    m = sa.atom_arrays.mass[:]
    v = sa.atom_arrays.v[:]
    n = length(sa.atom_arrays.r[:])
    vec = mean(sa.atom_arrays.v[:])
    @inbounds for i = 1:3
        @inbounds for j = 1:3
            @inbounds for k = 1:n
                kinetic_energy[i, j] +=
                    m[k] * (v[k][i] - vec[i]) * (v[k][j] - vec[j])
            end
        end
    end
    #println("kinetic energy $(kinetic_energy)")
    return 0.5 * kinetic_energy
end

# compute pressure using the atomic virial tensor and kinetic energy
function compute_pressure(
    simulation_arrays::SimulationArrays,
    cutoff,
    box_size,
    ::atomic_ke,
)
    #https://manual.gromacs.org/documentation/2021/reference-manual/algorithms/molecular-dynamics.html
    #https://manual.gromacs.org/documentation/2021/reference-manual/functions/long-range-vdw.html#virial
    #https://en.wikipedia.org/wiki/Virial_stress
    #https://aip.scitation.org/doi/pdf/10.1063/1.2363381
    vdw_virial_tensor = lj_virial_tensor(simulation_arrays, cutoff, box_size)
    bond_tensor = total_virial_bond_tensor(simulation_arrays, box_size)
    ke_tensor = kinetic_energy_tensor(simulation_arrays)
    vol = volume(box_size)
    total_tensor = (2 * ke_tensor - (vdw_virial_tensor + bond_tensor)) / vol
    press = tr(total_tensor) / 3.0
    return press #, (vdw_virial_tensor + bond_tensor)

end

# calculates the pressure tensor from the virial tensor
function press_tensor(vir_tensor::MMatrix, temp, n, vol)
    press = -vir_tensor / vol
    @inbounds for i = 1:3
        press[i, i] = press[i, i] + (n * 0.008314461 * temp / vol)
    end
    return press
end

# compute the pressure using the atomic virial and ideal gas
function compute_pressure(
    simulation_arrays::SimulationArrays,
    cutoff,
    box_size,
    temp,
    ::atomic_ig,
)
    virial_tensor = lj_virial_tensor(simulation_arrays, cutoff, box_size)
    bond_tensor = total_virial_bond_tensor(simulation_arrays, box_size)
    virial_tensor += bond_tensor
    vol = volume(box_size)
    n = length(simulation_arrays.atom_arrays.r[:])
    press = tr(press_tensor(virial_tensor, temp, n, vol)) / 3.0
    return press
end

# compute the pressure using the molecular virial and ideal gas
function compute_pressure(
    simulation_arrays::SimulationArrays,
    cutoff,
    box_size,
    temp,
    ::molecular_ig,
)
    n = length(simulation_arrays.molecule_arrays.COM[:])
    vol = volume(box_size)
    lj_virial_scalar = lj_virial(simulation_arrays, cutoff, box_size)
    pressure_ideal_gas = n * 0.0083144621 / vol * temp
    pressure_residual = -lj_virial_scalar / vol / 3.0
    return pressure_ideal_gas + pressure_residual
end

volume(box_size::SVector) = box_size[1] * box_size[2] * box_size[3]
press_full(vir::Float64, n::Int64, vol::Float64, tmp::Float64) =
    n * 0.0083144621 / vol * tmp + vir / vol / 3
kinetic_energy(v::Vector, m::Vector) = 0.5 * dot(m .* v, v)
