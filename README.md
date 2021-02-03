# AutoDiff_md

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BradenDKelly.github.io/AutoDiff_md.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BradenDKelly.github.io/AutoDiff_md.jl/dev)
[![Build Status](https://github.com/BradenDKelly/AutoDiff_md.jl/workflows/CI/badge.svg)](https://github.com/BradenDKelly/AutoDiff_md.jl/actions)
[![Coverage](https://codecov.io/gh/BradenDKelly/AutoDiff_md.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BradenDKelly/AutoDiff_md.jl)

This is a playground for the many features of Julia. It is meant as a place to learn/implement package development and feature development, such as using autodiff.

An example for using it is:

```Julia
temp = 119.8 * 2     # temperature (used for initial velocity)
dt = 0.002      # time step dimensionless
nsteps = 10000   # number of steps to use for trajectory
natoms = 1000
cutoff = 1.2 # nm
verlet_buffer = 0.2 #nm

top_file = joinpath(
    dirname(pathof(AutoDiff_md)),
    "../",
    "structures",
    "topology_files",
    "Ar.top",
    #"mea_tip3p.top",
) # Ar.top
specie_location = joinpath(
    dirname(pathof(AutoDiff_md)),
    "../",
    "structures",
    "single_molecules",
)  # "mea.pdb",
system_location =
    joinpath(dirname(pathof(AutoDiff_md)), "../", "structures", "whole_systems")

#records names, identifiers, coords of all atoms in system
#@time systemPDB = ReadPDB(joinpath(system_location, "mea_water.pdb")) # Ar_1000.pdb
@time systemPDB = ReadPDB(joinpath(system_location, "Ar_1000.pdb"))
#records all FF parameters for the system
@time systemTop = ReadTopFile(top_file) # returns struct FFParameters # located in Setup.jl
#specieList = ["mea.pdb", "tip3p.pdb"] # "Ar.pdb",
specieList = ["Ar.pdb"] # "Ar.pdb",

molecule_list = []
# records all names, identifiers, coords for all starting body_fixed molecules
for mol in specieList
    push!(molecule_list, ReadPDB(joinpath(specie_location, mol)))
end

atom_arrays, molecule_arrays =
    MakeAtomArrays(systemTop, systemPDB, molecule_list)  # located in setup.jl
intraFF, vdwTable, nonbonded_matrix, scaled_pairs =
    MakeTables(systemTop, systemPDB, molecule_list) # located in Setup.jl
box_size = SVector(systemPDB.box...)
verlet_list = VerletList(verlet_buffer, [], [])
make_neighbor_list!(
    atom_arrays.r,
    nonbonded_matrix,
    box_size,
    cutoff,
    verlet_list,
)

# this stucture holds information on the number of atoms, molecules, atom types, molecule types, charges
numbers = Numbers(
    length(atom_arrays),
    length(molecule_arrays),
    length(systemTop.atomTypes),
    length(systemTop.molParams),
    count(!iszero, atom_arrays.qq),
)
simulation_arrays = SimulationArrays(
    atom_arrays,
    molecule_arrays,
    intraFF,
    vdwTable,
    nonbonded_matrix,
    scaled_pairs,
    numbers,
    verlet_list,
)

simulation_controls = SimulationControls(
    VelocityVerlet(dt),
    AndersenThermostat(1.0, temp),
    MonteCarloBarostat(1.0, 0, 0, 14.0, volume(box_size) * 0.01, true),
)
mutable struct Properties
    kinetic_temp::Real
    pressure::Real
    total_energy::Real
end
props = Properties(0.0, 0.0, 0.0)
# equilibrate
start = Dates.now()
Logo()
# @btime analytical_total_force(simulation_arrays, cutoff, box_size)
# @code_llvm analytical_total_force(simulation_arrays, cutoff, box_size)
# @code_warntype analytical_total_force(simulation_arrays, cutoff, box_size)
simulate!(
    simulation_arrays,
    simulation_controls,
    box_size,
    nsteps,
    cutoff,
    props,
)

finish = Dates.now()
difference = finish - start
println(
    "Total runtime was: ",
    Dates.canonicalize(Dates.CompoundPeriod(difference)),
)


Completion()
