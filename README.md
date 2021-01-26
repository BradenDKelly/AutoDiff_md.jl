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

temp = 119.8    # temperature, Kelvin (used for initial velocity)
dt = 0.005      # time step dimensionless
dt = 0.0002     # time step with units (picoseconds)
nsteps = 3000   # number of steps to use for trajectory
natoms = 256

# parameters for Argon
epsilon, sigma, mass = 0.997, 0.3405, 39.94
eps = [epsilon for i=1:natoms]
sig = [sigma for i=1:natoms]
m = [mass for i=1:natoms]

box = 3.0 # nm
rho = natoms / box^3
# initial coordinates, r. Atoms placed on a lattice to ensure clean start.
r, L = initCubicGrid(natoms, rho)
box_size = SVector{3}(L, L, L)
cutoff = box_size[1] / 2 * 0.9
v = [velocity(mass, temp) for i in 1:natoms]
println("Starting box size is: $box_size")

# equilibrate
simulate(r, v, m, eps, sig, box_size, temp, dt, nsteps, cutoff)

# Production run
dt = 0.00109
nsteps = 30000
simulate(r, v, m, eps, sig, box_size, temp, dt, nsteps, cutoff)
```
