export MonteCarloBarostat, NoBarostat, apply_barostat!

include("types.jl")

#"""Struct with Monte Carlo style barostat parameters"""
mutable struct MonteCarloBarostat{F,I} <: Barostat
    pcouple::F
    vol_attempt::I
    vol_accept::I
    prev_attempt::I
    prev_accept::I
    update::I
    acceptance_ratio::F
    set_press::F
    vmax::F
    use::Bool
end

# TODO change eps, sig to atom_arrays and vdwTable
#"""Apply Monte Carlo style barostat"""
function apply_barostat!(
    simulation_array::SimulationArrays,
    barostat::MonteCarloBarostat,
    box_size,
    cutoff,
    temp,
)
    r, com, box_size, cutoff, nl =
        volume_change_lj_MC!(simulation_array, box_size, cutoff, temp, barostat)
    if barostat.vol_attempt % barostat.update == 0
        optimize_volume_change_step!(barostat, volume(box_size)*0.1)
    end
    #println("$(barostat.vmax), $(barostat.vol_attempt)")
    return r, com, box_size, cutoff, nl
end

#"""Struct for no barostat"""
struct NoBarostat <: Barostat
    use::Bool
end

#"""Apply no barostat"""
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

# optimize the vmax (max volume change) parameter
function optimize_volume_change_step!(barostat, L)

    # This code adjusts the maximum translation such that the probability of a successful move equals move_accept (fraction)
    # Written by Braden Kelly June 19, 2016 / modified for Julia April 2020
    # Based on the code by Frenkel and Smith (Fortran 77/90)
    #=
    ---------------------------
        naccepp is number of moves accepted as of previous call
        naccept is number of moves accepted as of this call
        attempp is number of attempts as of previous call
        attemt is number of attempts as of this call
        ratio is a numerical derivative
    =#
    if (barostat.prev_accept == 0) #! .or. attempp .ge. attempt
        barostat.prev_accept = barostat.vol_accept
        barostat.prev_attempt = barostat.vol_attempt
    else
        # numerical derivative
        ratio =
            real(barostat.vol_accept - barostat.prev_accept) /
            real(barostat.vol_attempt - barostat.prev_attempt)
        dr_old = barostat.vmax
        barostat.vmax = barostat.vmax * ratio / barostat.acceptance_ratio #move_accept
        dr_ratio = barostat.vmax / dr_old
        #! Place limits on changes
        if (dr_ratio > 1.5)
            barostat.vmax = dr_old * 1.5
        end
        if (dr_ratio < 0.5)
            barostat.vmax = dr_old * 0.5
        end
        if (barostat.vmax > L / 2)
            barostat.vmax = L / 2
        end
        barostat.prev_accept = barostat.vol_accept
        barostat.prev_attempt = barostat.vol_attempt
    end
    # If you are seeing this, then perhaps it is too late.
end # volume change
