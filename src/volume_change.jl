export volume_change_lj_MC!

#=
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
=#
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

#""" returns max and min values in an array of SVectors"""
function maxmin(array::Vector)
    maxVal = maximum(array[1])
    minVal = minimum(array[1])

    for svector in array
        if maximum(svector) > maxVal
            maxVal = maximum(svector)
        end
        if minimum(svector) < minVal
            minVal = minimum(svector)
        end
    end

    return minVal, maxVal
end

#=
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
=#
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
    #rememb = deepcopy(simulation_arrays.atom_arrays.r[1])
    #rememb2 = deepcopy(simulation_arrays.atom_arrays.r[2])
    old_box = deepcopy(box_size)
    #L_old = box_size[1]
    #volume_old = box_size[1] * box_size[2] * box_size[3]

    # calculate energy of original system
    energy_old = total_energy(simulation_arrays, cutoff, box_size, point, list)[1]

    # calculate random volume change
    change_vol = (rand() - 0.5) * barostat.vmax # note vmax is actually ln(vmax)
    #ln_vol_new = volume_old + change_vol
    #volume_new = ln_vol_new
    #L_new = volume_new^(1 / 3)
    #box_size = SVector(L_new, L_new, L_new)
    # update centers-of-mass coords
    simulation_arrays.molecule_arrays.COM[:] = total_com_update!(simulation_arrays)
    old_com = deepcopy(simulation_arrays.molecule_arrays.COM)
    old_r = deepcopy(simulation_arrays.atom_arrays.r[:])
    # update center of mass of molecules
    #new_atoms_positions, new_com, box_size = change_volume_molecular(simulation_arrays, change_vol, box_size)
    new_atoms_positions, new_com, box_size = change_volume_atomic(simulation_arrays, change_vol, box_size)
    #new_com = deepcopy(old_com .* (L_new / L_old))
    simulation_arrays.molecule_arrays.COM[:] = new_com
    simulation_arrays.atom_arrays.r[:] = new_atoms_positions

    # update neighborlist
    make_neighbor_list!(
        simulation_arrays.atom_arrays.r[:],
        simulation_arrays.nonbonded_matrix,
        box_size,
        cutoff,
        simulation_arrays.neighborlist
        )
    point = simulation_arrays.neighborlist.point[:]
    list = simulation_arrays.neighborlist.list[:]
    energy_new = total_energy(simulation_arrays, cutoff, box_size, point, list)[1]
    du = (energy_new - energy_old)
    dv = (volume(box_size) - volume(old_box))
    @assert isapprox(dv, change_vol)
    # walking in V - change (n) to (n+1) if walking in lnV
    test =
        -1 / (0.00831446 * temp) * (
            du + barostat.set_press * dv -
            (n) * 0.00831446 * temp * log(volume(box_size) / volume(old_box))
        )
    # test if we keep this volume or the original volume
    if rand() > exp(test)
        # reset to old coords, move failed
        simulation_arrays.atom_arrays.r[:] = old_r[:]
        simulation_arrays.molecule_arrays.COM[:] = old_com[:]
        box_size = old_box
        make_neighbor_list!(
            simulation_arrays.atom_arrays.r[:],
            simulation_arrays.nonbonded_matrix,
            box_size,
            cutoff,
            simulation_arrays.neighborlist
            )
    else
        # move was accepted, keep new coordinates and volume
        barostat.vol_accept += 1
        simulation_arrays.molecule_arrays.COM[:] = new_com[:]
    end
    # return (does not modify in place for some reason.)

    return simulation_arrays.atom_arrays.r[:], simulation_arrays.molecule_arrays.COM[:], box_size, cutoff, simulation_arrays.neighborlist

    # TODO update tail corrections
    # TODO add module for optimizing vmax
end

# molecular scaling for volume change
function change_volume_molecular(simulation_arrays::SimulationArrays, change, box_size)
    # change should be change in volume

    if length(change) == 1
        new_length = (volume(box_size) + change)^(1/3)
        #change = sign(change) * abs(change)^(1/3)    # convert to change in length
        #change = SVector(change, change, change) # assume same change in x, y, z
    end

    old_com = simulation_arrays.molecule_arrays.COM
    new_box = SVector(new_length, new_length, new_length)#box_size .+ change
    scaling = new_box ./ box_size

    # update centers of mass
    new_com = [xyz .* scaling for xyz in old_com]
    com_change = [new_com[i] .- old_com[i] for i=1:length(new_com)]
    # println("new scaling is $(L_new / L_old)")
    # for i = 1:length(new_com)
    #     println("i $i newcom $(new_com[i]) oldcom $(old_com[i])")
    # end

    # update atom coordinates
    # old_bonds = []
    # for i=1:2:n
    #     dist = norm(simulation_arrays.atom_arrays.r[i]- simulation_arrays.atom_arrays.r[i+1])
    #     push!(old_bonds, dist)
    #     if dist > 0.2
    #         println("SHIT particles $i, $(i+1) too far apart, $dist")
    #     end
    # end
    # for troubleshooting checks that particles are in box
    minV, maxV = maxmin(new_com)
    #println(minV, maxV)
    if minV < 0.0
        println("Shit, particle is less than 0")
    end
    if maxV > new_box[1]
        println("Shit, particle is outside box")
    end
    atom_coords = [zeros(SVector{3}) for i=1:length(simulation_arrays.atom_arrays.r[:])]
    for this_mol = 1:length(new_com)
        ia = simulation_arrays.molecule_arrays[this_mol].firstAtom
        ja = simulation_arrays.molecule_arrays[this_mol].lastAtom
        # TODO should not need a loop here, a vectorized fuse should work.
        for j = ia:ja
            atom_coords[j] = simulation_arrays.atom_arrays.r[j] .+ com_change[this_mol]
        end
    end
    #println("$(simulation_arrays.atom_arrays.r[1]), old was $rememb")
    #println("$(simulation_arrays.atom_arrays.r[2]), old was $rememb2")
    #sleep(180)
    #new_bonds = []
    for i=1:2:length(atom_coords)
        dist = norm(atom_coords[i]- atom_coords[i+1])
        #push!(new_bonds, dist)
        if abs(dist) > 0.2
            println("SHIT particles $i, $(i+1) too far apart, $dist")
        end
    end
    return atom_coords, new_com, new_box
end

function change_volume_atomic(simulation_arrays::SimulationArrays, change, box_size)
    # change should be change in volume

    if length(change) == 1
        new_length = (volume(box_size) + change)^(1/3)
        #change = sign(change) * abs(change)^(1/3)    # convert to change in length
        #change = SVector(change, change, change) # assume same change in x, y, z
    end
    old_length = volume(box_size)^(1/3)
    scale = new_length / old_length
    atom_coords = simulation_arrays.atom_arrays.r[:] .* scale
    new_com = total_com_update!(simulation_arrays)
    new_box = box_size .* scale
    return atom_coords, new_com, new_box
end
