#using StaticArrays
export
    pb!,
    integrator!

function pb!(r, box)
    """ Periodic Boundary Conditions

    Parameters
    ----------
    r : SVector{3}
        coordinates of a single atom
    box : Float64
        box size

    Returns
    Nothing. Coordinates are changed in place
    """
    x = r[1]
    y = r[2]
    z = r[3]

    if x < 0 x += box end
    if x > box x -= box end
    if y < 0 y += box end
    if y > box y -= box end
    if z < 0 z += box end
    if z > box z -= box end

    r = SVector{3}(x, y, z)
end

function integrator!(r, v, f, m, dt, eps, sig, cutoff, box_size, temp, thermo_couple)
    """
    Velocity-Verlet integrator scheme.

    Parameters
    ----------
    r : Vector{SVector}
        atom coordinates
    v : Vector{SVector}
        atom velocities
    f : Vector{SVector}
        forces on atoms
    m : Vector{Float64}
        masses of atoms
    dt : Float64
        time step (appx 0.005 for reduced units)
    eps : Vector
        epsilon parameters of all atoms
    sigma : Vector
        sigma parameters of all atoms
    cutoff : Float
        interaction cutoff distance
    box_size : SVector{3}
        vector of x, y, z box lengths (should all be the same)

    Returns
    Nothing. Updated in place.
    """
    n = length(r)
    for i=1:n
        # update velocities
        v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
        # update positons (apply periodic boundary conditions)
        r[i] = pb!(r[i] .+ v[i] .* dt, box_size[1])
        if maximum(r[i]) > box_size[1]
            println("atom $i is at coord $(r[i])")
        end
    end
    # update forces
    f = analytical_total_force(r, eps, sig, cutoff, box_size)
    for i=1:n
        # update velocities
        if rand() < thermo_couple * dt
            # Andersen Thermostat
            v[i] = velocity(m[i], temp)
        else
            v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
        end # if statement
    end # for loop
end # function
