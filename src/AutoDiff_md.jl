using BenchmarkTools
using LinearAlgebra
using StaticArrays
using Setfield
using Statistics

#include("energies.jl")
#include("forces.jl")
include("io_data.jl")
include("simulation.jl")
include("create_system.jl")



println("##########################################")
println("#                                        #")
println("#             New Simulation             #")
println("#                                        #")
println("##########################################")


# This simulation assumes mass = 1
temp = 1.       # temperature (not actually used since this is NVE)
dt = 0.005      # time step
nsteps = 10000  # number of steps to use for trajectory
# get initial coordinates, velocity, boxsize from file
r, v, L = ReadCNF("cnf.inp")
box_size = SVector{3}(L, L, L)
println("Starting box size is: $box_size")

simulate(r, v, box_size, temp, dt, nsteps)

println("##########################################")
println("#                                        #")
println("#        Simulation Complete             #")
println("#                                        #")
println("##########################################")


#module AutoDiff_md

# Write your package code here.

#end


# number of atoms
#natoms = 32    # if using cnf this is over-ridden and is 256
#startingConfiguration = "cnf"

# coordinates of the atoms. r is a list, each element is a static vector len=3
#if startingConfiguration == "random"
#    # Note, we will use square box, this is just a formality
#    b = 17.0
#    box_size = SVector{3}(b, b, b)
#    r = [SVector{3}(rand(), rand(), rand()) .* box_size[1] for i = 1:natoms]
#    println("ρ = $(natoms / box_size[1]^3)")
#elseif startingConfiguration == "cnf"
#    r, e, L = ReadCNF("cnf.inp")
#    box_size = SVector{3}(L, L, L)
#    println("Starting box size is: $box_size")
#else
#    rho = 0.01
#    r, L = InitCubicGrid(natoms, rho)
#    box_size = SVector{3}(L, L, L)
#    println("Starting box size is: $box_size")
#end
