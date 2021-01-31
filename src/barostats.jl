export
    MonteCarloBarostat,
    NoBarostat,
    apply_barostat!

mutable struct MonteCarloBarostat{F, I} <: Barostat
    pcouple::F
    vmax::F
    vol_attempt::I
    vol_accept::I
end

# TODO change eps, sig to atom_arrays and vdwTable
function apply_barostat!(r, box_size, cutoff, eps, sig, temp, press, vmax, barostat::MonteCarloBarostat)
    r, box_size, cutoff = volume_change_lj_MC!(r, box_size, cutoff, eps, sig, temp, press, vmax)
end

struct NoBarostat end

function apply_barostat!(::NoBarostat)
    return None
end
