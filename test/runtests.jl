using AutoDiff_md
using LinearAlgebra
using Test
using StaticArrays
using BenchmarkTools

LJ_energy(eps, sig, r) = 4 * eps *( (sig / r)^12 - (sig / r)^6)
LJ_force(eps, sig, r) = 48 * eps * r / r^2 * ( (sig / r)^12 - 0.5 * (sig / r)^6)
@testset "AutoDiff_md.jl" begin
    # get initial coordinates, velocity, boxsize from file
    r, v, L = ReadCNF("cnf.inp")
    @test length(r) == 256
    @test length(r) != 255
    @test length(v) == 256
    @test length(v) != 257
    @test typeof(L) == Float64
    @test typeof(L) != Int64
    @test typeof(r) ==Array{SVector{3, Float64}, 1}

    box = SVector(3.0, 3.0, 3.0)
    cut = 1.4
    r1 = SVector(0.0, 0.0, 0.0)
    r2 = SVector(1.0, 0.0, 0.0)
    eps, sig = 1.0, 1.0
    # test pair energy
    # @test LJ_energy(eps, sig, 1.0) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # # next should fail, rij longer than cutoff
    # r2 = SVector(2.0, 0.0, 0.0)
    # @test LJ_energy(eps, sig, 2.0) != pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # # next should pass, just inside cutoff
    # r2 = SVector(0.0, 1.39, 0.0)
    # @test LJ_energy(eps, sig, 1.39) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # # next should fail, just outside cutoff
    # r2 = SVector(0.0, 0.0, 1.41)
    # @test LJ_energy(eps, sig, 1.41) != pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    #
    # # test pair_force
    # r2 = SVector(1.0, 0.0, 0.0)
    # @test LJ_force(eps, sig, 1.0) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # # next should fail, rij longer than cutoff
    # r2 = SVector(2.0, 0.0, 0.0)
    # @test LJ_force(eps, sig, 2.0) !=  grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # # next should pass, just inside cutoff
    # r2 = SVector(0.0, 1.39, 0.0)
    # @test LJ_force(eps, sig, 1.39) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[2]
    # next should fail, just outside cutoff
    r2 = SVector(0.0, 0.0, 1.41)
    @test LJ_force(eps, sig, 1.41) != grad(r2, r1, eps, eps, sig, sig, cut, box)[3]
    @test LJ_force(eps, sig, 1.41) != lj_grad(r2, r1, eps, sig, cut, box)[3]

    #
    # eps, sig = 0.95, 0.34
    # # test pair energy
    # r2 = SVector(1.0, 0.0, 0.0)
    # @test LJ_energy(eps, sig, 1.0) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # # next should fail, rij longer than cutoff, but not mirror image
    # r2 = SVector(2.0, 0.0, 0.0)
    # @test LJ_energy(eps, sig, 2.0) !=  pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # # next should pass, shorter than cutoff, shorter than mirror image
    # r2 = SVector(1.0, 0.0, 0.0)
    # @test LJ_force(eps, sig, 1.0) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # # next should fail, rij longer than cutoff, shorter than mirror image
    # r2 = SVector(0.0, 1.45, 0.0)
    # @test LJ_force(eps, sig, 1.45) !=  grad(r2, r1, eps, eps, sig, sig, cut, box)[2]

    # test total_energy
    """
    The following hardcoded values were calculated using the Python code
    that is companion to Allen & Tildesley Computer Simulation of Liquids, 2017

    Their Python code is numpy accelerated and super slow. Don't vote Python.

    Also, fuck Python's 0-indexing.
    """
    # natoms = length(r)
    # box_size = SVector(L, L, L)
    # cut = 2.5
    # eps, sig = 1, 1
    # eps = [eps for i=1:natoms]
    # sig = [sig for i=1:natoms]
    # m = [1.0 for i=1:natoms]
    # # make momentum center of mass equal to 0 i.e., our box is not moving
    # v = initial_velocity(m, v)
    # # test kinetic energy
    # KE = kinetic_energy(v, m)
    # @test isapprox(KE, 85.33333333639997, atol=1e-6, rtol=1e-6)
    # temp = 2 * kinetic_energy(v, m) / (3 * natoms - 3)
    # @test  isapprox(temp, 0.223094, atol=1e-6, rtol=1e-6)
    # # test potential energy
    # PE = total_energy(r, eps, sig, cut, box_size)[1]
    # @test isapprox(PE, -318.00164297886465, atol=1e-6, rtol=1e-6)
    # # test total energy
    # @test isapprox(KE + PE, -232.668416, atol=1e-6, rtol=1e-6)
    # # test total_force
    # f = analytical_total_force(r, eps, sig, cut, box_size)
    # fsq = sum(sum(x->x.^2, f))
    # @test isapprox(fsq, 3.1205876093428077e-14, atol=1e-18, rtol=1e-18)

    # f1 = SVector(1.27370661e-08, 1.27370659e-08, 1.27370659e-08)
    # f100 = SVector(-3.33066907e-16,  3.33066907e-16,  5.45664347e-10)
    # f200 = SVector(-1.27370662e-08,  1.27370661e-08,  0.00000000e+00)
    # f256 = SVector(-1.27370662e-08, -5.45664347e-10, -1.27370661e-08)
    # @test isapprox(f[1], f1, atol=1e-12, rtol=1e-12)
    # @test isapprox(f[100], f100, atol=1e-14, rtol=1e-14)
    # @test isapprox(f[200], f200, atol=1e-12, rtol=1e-12)
    # @test isapprox(f[256], f256, atol=1e-12, rtol=1e-12)
    # @test !isapprox(f[256], f256, atol=1e-16, rtol=1e-16)
    # # test 4, 32, 64, 216, 512, 1024 particle size
    #
    # # test periodic boundary conditions. assumes box is 0 -> box_size
    # @test isapprox(pb!(SVector(2.0, 0.0, 0.0), 1.5), SVector(0.5, 0.0, 0.0))
    # @test !isapprox(pb!(SVector(2.0, 0.0, 0.0), 1.5), SVector(2.0, 0.0, 0.0))
    # @test isapprox(pb!(SVector(-0.5, 0.0, 0.0), 1.5), SVector(1.0, 0.0, 0.0))
    #
    # @test isapprox(pb!(SVector(0.0, 2.0, 0.0), 1.5), SVector(0.0, 0.5, 0.0))
    # @test isapprox(pb!(SVector(0.0, -1.0, 0.0), 1.5), SVector(0.0, 0.5, 0.0))
    #
    # @test isapprox(pb!(SVector(0.0, 0.0, 2.0), 1.5), SVector(0.0, 0.0, 0.5))
    # @test isapprox(pb!(SVector(0.0, 0.0, -1.3), 1.5), SVector(0.0, 0.0, 0.2))
    # test maxwellboltzmann

    # test input read

    # test integrators

    # test thermostats

    # test barostats

    # test energy conservation
        # run short run, make sure std dev is small...?

    # test pressure

    # temp = 1.0     # temperature (used for initial velocity)
    # dt = 0.005      # time step dimensionless
    # nsteps = 7000
    # mutable struct Properties
    #     kinetic_temp::Real
    #     pressure::Real
    #     total_energy::Real
    # end
    # props = Properties(0.0, 0.0, 0.0)
    # # TODO add @test_logs or suppressor.jl here to avoid print output
    # simulate!(r, v, m, eps, sig, box_size, temp, dt, nsteps, cut, props)
    # @test  0.90 < props.kinetic_temp * 0.00831446 < 1.1
    # @test  -0.1 < props.pressure < 0.0
    # @test  -250 < props.total_energy < -240

    ######### DIATOMICS ###########
    top_file = joinpath(
        dirname(pathof(AutoDiff_md)),
        "../",
        "structures",
        "topology_files",
        "N2_lj.top"
    )
    specie_location = joinpath(
        dirname(pathof(AutoDiff_md)),
        "../",
        "structures",
        "single_molecules",
    )
    system_location =
        joinpath(dirname(pathof(AutoDiff_md)), "../", "structures", "whole_systems")
    systemPDB = ReadPDB(joinpath(system_location, "N2_lj_1000.pdb"))
    #records all FF parameters for the system
    systemTop = ReadTopFile(top_file)
    @test systemTop.defaults.nbfunc == 1
    @test systemTop.defaults.combRule == 2
    @test systemTop.defaults.genPairs == "no"
    @test systemTop.defaults.fudgeLJ == 1.00
    @test systemTop.defaults.fudgeQQ == 1.00
    @test systemTop.atomTypes[1].name == "A1"
    @test systemTop.atomTypes[2].name == "A2"
    @test systemTop.atomTypes[1].atomicnr == "A1"
    @test systemTop.atomTypes[1].mass == 14.007
    @test systemTop.atomTypes[2].mass == 14.007
    @test systemTop.atomTypes[1].charge == 0.00
    @test systemTop.atomTypes[1].σ == 0.3262
    @test systemTop.atomTypes[1].ϵ == 0.3177
    @test systemTop.molParams[1].name == "N2"
    @test systemTop.molParams[1].nrexcl == 1
    @test length(systemTop.molParams[1].atoms) == 2
    @test systemTop.molParams[1].atoms[1].type == "A1"
    @test systemTop.molParams[1].atoms[1].resnr == 1
    @test systemTop.molParams[1].atoms[1].nr == 1
    @test systemTop.molParams[1].atoms[1].resnm == "N2"
    @test systemTop.molParams[1].atoms[1].atomnm == "A1"
    @test systemTop.molParams[1].atoms[1].mass == 14.007
    @test systemTop.molParams[1].atoms[1].charge == 0.00
    @test length(systemTop.molParams[1].bonds) == 1
    @test systemTop.molParams[1].bonds[1].ai == 1
    @test systemTop.molParams[1].bonds[1].aj == 2
    @test systemTop.molParams[1].bonds[1].funct == 1
    @test systemTop.molParams[1].bonds[1].bondLength == 0.147
    @test systemTop.molParams[1].bonds[1].bondLength != 0.149
    @test systemTop.molParams[1].bonds[1].kparam == 268280.0
    @test systemTop.molParams[1].bonds[1].kparam != 2268280.0
    @test length(systemTop.molParams[1].angles) == 0
    @test length(systemTop.molParams[1].angles) != 1
    @test length(systemTop.molParams[1].dihedrals) == 0
    @test length(systemTop.molParams[1].dihedrals) != 1

    for mol in specieList
        push!(molecule_list, ReadPDB(joinpath(specie_location, mol)))
    end
    atom_arrays, molecule_arrays =
        MakeAtomArrays(systemTop, systemPDB, molecule_list)
    intraFF, vdwTable, nonbonded_matrix, scaled_pairs =
        MakeTables(systemTop, systemPDB, molecule_list) # located in Setup.jl
    box_size = SVector(systemPDB.box...)
    verlet_list = VerletList(verlet_buffer, [], [], 15)
    make_neighbor_list!(
        atom_arrays.r,
        nonbonded_matrix,
        box_size,
        cutoff,
        verlet_list,
    )
    simulation_arrays = SimulationArrays(
        atom_arrays,
        molecule_arrays,
        verlet_list,
        intraFF,
        vdwTable,
        nonbonded_matrix,
        scaled_pairs,
        numbers,
    )
    @test length(simulation_arrays.atom_arrays) == 2000
    @test length(simulation_arrays.molecule_arrays) == 1000
    @test size(simulation_arrays.nonbonded_matrix) == (2000,2000)

    @test simulation_arrays.atom_arrays[1].molNum == 1
    @test simulation_arrays.atom_arrays[2].molNum == 1
    @test simulation_arrays.atom_arrays[3].molNum == 2
    @test simulation_arrays.atom_arrays[1998].molNum == 999
    @test simulation_arrays.atom_arrays[1999].molNum == 1000
    @test simulation_arrays.atom_arrays[2000].molNum == 1000
    @test simulation_arrays.atom_arrays[2].molNum != 2

    @test simulation_arrays.atom_arrays[1].molType == 1
    @test simulation_arrays.atom_arrays[2].molType == 1
    @test simulation_arrays.atom_arrays[3].molType == 1
    @test simulation_arrays.atom_arrays[1998].molType == 1
    @test simulation_arrays.atom_arrays[1999].molType == 1
    @test simulation_arrays.atom_arrays[2000].molType == 1
    @test simulation_arrays.atom_arrays[2].molType != 2

    @test simulation_arrays.atom_arrays[1].atype == 1
    @test simulation_arrays.atom_arrays[2].atype == 2
    @test simulation_arrays.atom_arrays[3].atype == 1
    @test simulation_arrays.atom_arrays[1998].atype == 2
    @test simulation_arrays.atom_arrays[1999].atype == 1
    @test simulation_arrays.atom_arrays[2000].atype == 2
    @test simulation_arrays.atom_arrays[2].atype != 1

    @test simulation_arrays.atom_arrays[1].mass == 14.007
    @test simulation_arrays.atom_arrays[2].mass == 14.007
    @test simulation_arrays.atom_arrays[3].mass == 14.007
    @test simulation_arrays.atom_arrays[1998].mass == 14.007
    @test simulation_arrays.atom_arrays[1999].mass == 14.007
    @test simulation_arrays.atom_arrays[2000].mass == 14.007
    @test simulation_arrays.atom_arrays[2].mass != 12.007

    @test simulation_arrays.molecule_arrays[1].firstAtom == 1
    @test simulation_arrays.molecule_arrays[1].lastAtom == 2
    @test simulation_arrays.molecule_arrays[2].firstAtom == 3
    @test simulation_arrays.molecule_arrays[2].lastAtom == 4
    @test simulation_arrays.molecule_arrays[999].firstAtom == 1997
    @test simulation_arrays.molecule_arrays[999].lastAtom == 1998
    @test simulation_arrays.molecule_arrays[1000].firstAtom == 1999
    @test simulation_arrays.molecule_arrays[1000].lastAtom == 2000
    @test simulation_arrays.molecule_arrays[1].molType == 1
    @test simulation_arrays.molecule_arrays[1000].molType == 1
    @test simulation_arrays.molecule_arrays[1].mass == 28.014
    @test simulation_arrays.molecule_arrays[1000].mass == 28.014

    #= TODO
    test FindMolType and friends
    =#
    function set_norm(norm)
        # assumes one coord is at 0, 0, 0
        # i.e., sqrt ( x^2 + y^2 + z^2) = norm, but, x = y = z
        # sqrt(3*x^2) = norm ---> x = sqrt(norm^2 / 3 ) and x = y = z
        r = sqrt(norm^2 / 3)
        return SVector(r, r, r)
    end
    # Test bonds
    # choose coords such that norm of difference is 1
    coord1 = SVector(0.,0.,0.)
    coord2 = set_norm(1.0)
    box_size = SVector(3.,3.,3.)
    kspring = 1.0
    bond_length = 1.0
    @test bond_force(coord1, coord2, kspring, bond_length, box_size) == (SVector(0,0,0), SVector(0,0,0))

    # choose coords such that norm of difference is 0.5
    normal = 0.5
    coord2 = set_norm(normal)
    @test norm(vector(coord1, coord2, box_size)) == normal
    rab = vector(coord1, coord2, box_size)
    @test norm(normalize(rab)) == 1.0
    @test bond_force(coord1, coord2, kspring, bond_length, box_size) == (SVector(-normal * normalize(rab)... ), SVector(normal * normalize(rab)...))
    kspring = 2000.0
    @test bond_force(coord1, coord2, kspring, bond_length, box_size) == (SVector(-normal * kspring * normalize(rab)... ), SVector(normal * kspring * normalize(rab)...))
    # force mirror image separation to kick in
    box_size = SVector(0.3, 0.3, 0.3)
    @test bond_force(coord1, coord2, kspring, bond_length, box_size) != (SVector(-normal * kspring * normalize(rab)... ), SVector(normal * kspring * normalize(rab)...))
end
