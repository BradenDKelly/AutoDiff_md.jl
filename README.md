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

# This simulation assumes mass = 1
temp = 1.       # temperature (not actually used since this is NVE)
dt = 0.005      # time step
nsteps = 10000  # number of steps to use for trajectory

# get initial coordinates, velocity, boxsize from file
r, v, L = ReadCNF("cnf.inp")

box_size = SVector{3}(L, L, L)
m = ones(length(r)) # all masses are unit 1

# run the simulation
simulate(r, v, m, box_size, temp, dt, nsteps)
```
