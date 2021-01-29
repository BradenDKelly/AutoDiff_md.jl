# AutoDiff_md

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BradenDKelly.github.io/AutoDiff_md.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BradenDKelly.github.io/AutoDiff_md.jl/dev)
[![Build Status](https://github.com/BradenDKelly/AutoDiff_md.jl/workflows/CI/badge.svg)](https://github.com/BradenDKelly/AutoDiff_md.jl/actions)
[![Coverage](https://codecov.io/gh/BradenDKelly/AutoDiff_md.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BradenDKelly/AutoDiff_md.jl)

This is a playground for the many features of Julia. It is meant as a place to learn/implement package development and feature development, such as using autodiff.

An example for using it is:

```Julia
using AutoDiff_md
using StaticArrays

temp = 119.8 * 2     # temperature (used for initial velocity)
dt = 0.002      # time step dimensionless
nsteps = 50000   # number of steps to use for trajectory
natoms = 1000

# parameters for Argon
epsilon, sigma, mass = 0.996066, 0.3405, 39.96
eps = [epsilon for i=1:natoms]
sig = [sigma for i=1:natoms]
m = [mass for i=1:natoms]
r, v, L = Read_gromacs("nvt_equil.g96")
press = 14.0
pcouple = 10 # # time steps between vol changes

box_size = SVector{3}(L, L, L)
cutoff = 1.2 # box_size[1] / 2 * 0.9
v = [velocity(mass, temp) for i in 1:natoms]

mutable struct Properties
    kinetic_temp::Real
    pressure::Real
    total_energy::Real
end
props = Properties(0.0, 0.0, 0.0)
simulate!(r, v, m, eps, sig, box_size, temp, press, dt, nsteps, cutoff, props, pcouple)
