export
    MonteCarloBarostat,
    NoBarostat,
    apply_barostat!

include("types.jl")
mutable struct MonteCarloBarostat{F, I} <: Barostat
    pcouple::F
    vol_attempt::I
    vol_accept::I
    set_press::F
    vmax::F
    use::Bool
end

# TODO change eps, sig to atom_arrays and vdwTable
function apply_barostat!(
    simulation_array::SimulationArrays,
    barostat::MonteCarloBarostat,
    box_size,
    cutoff,
    temp
    )
    r, box_size, cutoff = volume_change_lj_MC!(simulation_array, box_size, cutoff, temp, barostat)
    return r, box_size, cutoff
end

struct NoBarostat <: Barostat
    use::Bool
end

function apply_barostat!(r, box_size, cutoff, eps, sig, temp, press, vmax, barostat::NoBarostat)
    return r, box_size, cutoff
end
