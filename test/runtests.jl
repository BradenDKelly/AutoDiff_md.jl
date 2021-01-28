using AutoDiff_md
using Test
using StaticArrays

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
    @test LJ_energy(eps, sig, 1.0) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # next should fail, rij longer than cutoff
    r2 = SVector(2.0, 0.0, 0.0)
    @test LJ_energy(eps, sig, 2.0) != pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # next should pass, just inside cutoff
    r2 = SVector(0.0, 1.39, 0.0)
    @test LJ_energy(eps, sig, 1.39) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # next should fail, just outside cutoff
    r2 = SVector(0.0, 0.0, 1.41)
    @test LJ_energy(eps, sig, 1.41) != pair_energy(r1, r2, eps, eps, sig, sig, cut, box)

    # test pair_force
    r2 = SVector(1.0, 0.0, 0.0)
    @test LJ_force(eps, sig, 1.0) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # next should fail, rij longer than cutoff
    r2 = SVector(2.0, 0.0, 0.0)
    @test LJ_force(eps, sig, 2.0) !=  grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # next should pass, just inside cutoff
    r2 = SVector(0.0, 1.39, 0.0)
    @test LJ_force(eps, sig, 1.39) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[2]
    # next should fail, just outside cutoff
    r2 = SVector(0.0, 0.0, 1.41)
    @test LJ_force(eps, sig, 1.41) != grad(r2, r1, eps, eps, sig, sig, cut, box)[3]

    eps, sig = 0.95, 0.34
    # test pair energy
    r2 = SVector(1.0, 0.0, 0.0)
    @test LJ_energy(eps, sig, 1.0) ≈ pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # next should fail, rij longer than cutoff, but not mirror image
    r2 = SVector(2.0, 0.0, 0.0)
    @test LJ_energy(eps, sig, 2.0) !=  pair_energy(r1, r2, eps, eps, sig, sig, cut, box)
    # next should pass, shorter than cutoff, shorter than mirror image
    r2 = SVector(1.0, 0.0, 0.0)
    @test LJ_force(eps, sig, 1.0) ≈ grad(r2, r1, eps, eps, sig, sig, cut, box)[1]
    # next should fail, rij longer than cutoff, shorter than mirror image
    r2 = SVector(0.0, 1.45, 0.0)
    @test LJ_force(eps, sig, 1.45) !=  grad(r2, r1, eps, eps, sig, sig, cut, box)[2]

    # test total_energy
    """
    The following hardcoded values were calculated using the Python code
    that is companion to Allen & Tildesley Computer Simulation of Liquids, 2017

    Their Python code is numpy accelerated and super slow. Don't vote Python.

    Also, fuck Python's 0-indexing.
    """
    natoms = length(r)
    box_size = SVector(L, L, L)
    cut = 2.5
    eps, sig = 1, 1
    eps = [eps for i=1:natoms]
    sig = [sig for i=1:natoms]
    m = [1.0 for i=1:natoms]
    # make momentum center of mass equal to 0 i.e., our box is not moving
    v = initial_velocity(m, v)
    # test kinetic energy
    KE = kinetic_energy(v, m)
    @test isapprox(KE, 85.33333333639997, atol=1e-6, rtol=1e-6)
    temp = 2 * kinetic_energy(v, m) / (3 * natoms - 3)
    @test  isapprox(temp, 0.223094, atol=1e-6, rtol=1e-6)
    # test potential energy
    PE = total_energy(r, eps, sig, cut, box_size)[1]
    @test isapprox(PE, -318.00164297886465, atol=1e-6, rtol=1e-6)
    # test total energy
    @test isapprox(KE + PE, -232.668416, atol=1e-6, rtol=1e-6)
    # test total_force
    f = analytical_total_force(r, eps, sig, cut, box_size)
    fsq = sum(sum(x->x.^2, f))
    @test isapprox(fsq, 3.1205876093428077e-14, atol=1e-18, rtol=1e-18)

    f1 = SVector(1.27370661e-08, 1.27370659e-08, 1.27370659e-08)
    f100 = SVector(-3.33066907e-16,  3.33066907e-16,  5.45664347e-10)
    f200 = SVector(-1.27370662e-08,  1.27370661e-08,  0.00000000e+00)
    f256 = SVector(-1.27370662e-08, -5.45664347e-10, -1.27370661e-08)
    @test isapprox(f[1], f1, atol=1e-12, rtol=1e-12)
    @test isapprox(f[100], f100, atol=1e-14, rtol=1e-14)
    @test isapprox(f[200], f200, atol=1e-12, rtol=1e-12)
    @test isapprox(f[256], f256, atol=1e-12, rtol=1e-12)
    @test !isapprox(f[256], f256, atol=1e-16, rtol=1e-16)
    # test 4, 32, 64, 216, 512, 1024 particle size

    # test periodic boundary conditions. assumes box is 0 -> box_size
    @test isapprox(pb!(SVector(2.0, 0.0, 0.0), 1.5), SVector(0.5, 0.0, 0.0))
    @test !isapprox(pb!(SVector(2.0, 0.0, 0.0), 1.5), SVector(2.0, 0.0, 0.0))
    @test isapprox(pb!(SVector(-0.5, 0.0, 0.0), 1.5), SVector(1.0, 0.0, 0.0))

    @test isapprox(pb!(SVector(0.0, 2.0, 0.0), 1.5), SVector(0.0, 0.5, 0.0))
    @test isapprox(pb!(SVector(0.0, -1.0, 0.0), 1.5), SVector(0.0, 0.5, 0.0))

    @test isapprox(pb!(SVector(0.0, 0.0, 2.0), 1.5), SVector(0.0, 0.0, 0.5))
    @test isapprox(pb!(SVector(0.0, 0.0, -1.3), 1.5), SVector(0.0, 0.0, 0.2))
    # test maxwellboltzmann

    # test input read

    # test integrators

    # test thermostats

    # test barostats

    # test energy conservation
        # run short run, make sure std dev is small...?

    # test pressure

    temp = 1.0     # temperature (used for initial velocity)
    dt = 0.005      # time step dimensionless
    nsteps = 7000
    mutable struct Properties
        kinetic_temp::Real
        pressure::Real
        total_energy::Real
    end
    props = Properties(0.0, 0.0, 0.0)
    # TODO add @test_logs or suppressor.jl here to avoid print output
    simulate!(r, v, m, eps, sig, box_size, temp, dt, nsteps, cut, props)
    @test  0.90 < props.kinetic_temp * 0.00831446 < 1.1
    @test  -0.1 < props.pressure < 0.0
    @test  -250 < props.total_energy < -240





end
