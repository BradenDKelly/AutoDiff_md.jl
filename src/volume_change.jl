export
    volume_change_lj_MC!

function volume_change_lj_MC!(r, box_size, eps, sig, cutoff, temp, press, vmax)
    """Changes the volume of the simulation box using stat mech"""
    energy_old = total_energy(r./box_size[1], eps, sig, cutoff/box_size[1], box_size./box_size[1])[1]

    L_old = box_size[1]
    volume_old = box_size[1] * box_size[2] * box_size[3]
    change_vol = (rand() - 0.5) * vmax # note vmax is actually ln(vmax)
    ln_vol_new = log(volume_old) + change_vol
    volume_new = exp(ln_vol_new)
    L_new = volume_new^(1/3)
    box_size = SVector(L_new, L_new, L_new)
    cutoff = min(cutoff, L_new / 2)
    n = length(r)
    r = r .* (L_new / L_old)

    energy_new = total_energy(r./box_size[1], eps, sig, cutoff/box_size[1], box_size./box_size[1])[1]
    du = (energy_new - energy_old)
    dv = (volume_new - volume_old)
    #@assert dv == exp(change_vol)

    test = -1 / (0.00831446 * temp) * ( du +
            press * dv -
            (n+1) * 0.00831446 * temp * log(volume_new / volume_old))
    #println("du $du, dv $dv, test $test")
    if rand() > exp(test)
        r = r .* (L_old / L_new)
        box_size = SVector(L_old, L_old, L_old)
    end
    return r, box_size
end
