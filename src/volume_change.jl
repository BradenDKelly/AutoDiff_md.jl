export volume_change_lj_MC!

"""Changes the volume of the simulation box using stat mech
Parameters
----------
r : Vector{SVector}
    atom coordinates
v : Vector{SVector}
    atom velocities
box_size : SVector{3}
    vector of x, y, z box lengths (should all be the same)
cutoff : Float
    interaction cutoff distance
eps : Vector
    epsilon parameters of all atoms
sigma : Vector
    sigma parameters of all atoms
temp : Float
    temperature of the simulation
press : Float
    pressure of the simulation
vmax : Float
    maximum volume change (max logarithm of volume)

Returns
r, boxsize, cutoff for the new volume
"""
function volume_change_lj_MC!(r, box_size, cutoff, eps, sig, temp, press, vmax)

    n = length(r)
    L_old = box_size[1]
    volume_old = box_size[1] * box_size[2] * box_size[3]

    # calculate energy of original system
    energy_old = total_energy(r, eps, sig, cutoff, box_size)[1]

    # calculate random volume change
    change_vol = (rand() - 0.5) * vmax # note vmax is actually ln(vmax)
    ln_vol_new = log(volume_old) + change_vol
    volume_new = exp(ln_vol_new)
    L_new = volume_new^(1 / 3)
    box_size = SVector(L_new, L_new, L_new)
    cutoff = min(cutoff, L_new / 2)
    r = r .* (L_new / L_old)
    # calculate energy of system with new volume
    energy_new = total_energy(r, eps, sig, cutoff, box_size)[1]

    du = (energy_new - energy_old)
    dv = (volume_new - volume_old)
    test =
        -1 / (0.00831446 * temp) * (
            du + press * dv -
            (n + 1) * 0.00831446 * temp * log(volume_new / volume_old)
        )

    # test if we keep this volume or the original volume
    if rand() > exp(test)
        r = r .* (L_old / L_new)
        box_size = SVector(L_old, L_old, L_old)
    end
    # return (does not modify in place for some reason.)
    return r, box_size, cutoff

    # TODO update tail corrections
end

"""Changes the volume of the simulation box using stat mech
    Parameters
    ----------
    simulation_array : SimulationArrays
        holds r, v, f for atoms, r for molecules and other parameters.
    box_size : SVector{3}
        vector of x, y, z box lengths (should all be the same)
    cutoff : Float
        interaction cutoff distance
    temp : Float
        temperature
    barostat : Barostat
        pcouple::F
        vol_attempt::I
        vol_accept::I
        set_press::F
        vmax::F
        use::Bool
    Returns
    r, boxsize, cutoff for the new volume
"""
function volume_change_lj_MC!(
    simulation_arrays::SimulationArrays,
    box_size,
    cutoff,
    temp,
    barostat::MonteCarloBarostat,
)
    point = simulation_arrays.neighborlist.point[:]
    list = simulation_arrays.neighborlist.list[:]

    barostat.vol_attempt += 1
    n = length(simulation_arrays.atom_arrays.r[:])
    L_old = box_size[1]
    volume_old = box_size[1] * box_size[2] * box_size[3]

    # calculate energy of original system
    energy_old = total_energy(simulation_arrays, cutoff, box_size, point, list)[1]

    # calculate random volume change
    change_vol = (rand() - 0.5) * barostat.vmax # note vmax is actually ln(vmax)
    ln_vol_new = log(volume_old) + change_vol
    volume_new = exp(ln_vol_new)
    L_new = volume_new^(1 / 3)
    box_size = SVector(L_new, L_new, L_new)
    cutoff = min(cutoff, L_new / 2)
    simulation_arrays.atom_arrays.r[:] =
        simulation_arrays.atom_arrays.r[:] .* (L_new / L_old)
    # calculate energy of system with new volume
    energy_new = total_energy(simulation_arrays, cutoff, box_size, point, list)[1]

    du = (energy_new - energy_old)
    dv = (volume_new - volume_old)
    test =
        -1 / (0.00831446 * temp) * (
            du + barostat.set_press * dv -
            (n + 1) * 0.00831446 * temp * log(volume_new / volume_old)
        )

    # test if we keep this volume or the original volume
    if rand() > exp(test)
        simulation_arrays.atom_arrays.r[:] =
            simulation_arrays.atom_arrays.r[:] .* (L_old / L_new)
        box_size = SVector(L_old, L_old, L_old)
    else
        barostat.vol_accept += 1
    end
    # return (does not modify in place for some reason.)
    return simulation_arrays.atom_arrays.r[:], box_size, cutoff

    # TODO update tail corrections
    # TODO add module for optimizing vmax
end
