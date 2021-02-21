export MonteCarloBarostat, NoBarostat, apply_barostat!

include("types.jl")

"""Struct with Monte Carlo style barostat parameters"""
mutable struct MonteCarloBarostat{F,I} <: Barostat
    pcouple::F
    vol_attempt::I
    vol_accept::I
    set_press::F
    vmax::F
    use::Bool
end

# TODO change eps, sig to atom_arrays and vdwTable
"""Apply Monte Carlo style barostat"""
function apply_barostat!(
    simulation_array::SimulationArrays,
    barostat::MonteCarloBarostat,
    box_size,
    cutoff,
    temp,
)
    r, com, box_size, cutoff =
        volume_change_lj_MC!(simulation_array, box_size, cutoff, temp, barostat)
    return r, com, box_size, cutoff
end

"""Struct for no barostat"""
struct NoBarostat <: Barostat
    use::Bool
end

"""Apply no barostat"""
function apply_barostat!(
    r,
    box_size,
    cutoff,
    eps,
    sig,
    temp,
    press,
    vmax,
    barostat::NoBarostat,
)
    # TODO figure out why I need to return, and update in place doesn't work
    return r, box_size, cutoff
end
