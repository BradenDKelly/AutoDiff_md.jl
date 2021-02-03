#using StaticArrays
export pb!, apply_integrator!, VelocityVerlet

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

    if x < 0
        x += box
    end
    if x > box
        x -= box
    end
    if y < 0
        y += box
    end
    if y > box
        y -= box
    end
    if z < 0
        z += box
    end
    if z > box
        z -= box
    end

    r = SVector{3}(x, y, z)
end

struct VelocityVerlet{F} <: Integrator
    dt::F
end

function apply_integrator!(
    simulation_arrays,
    integrator::VelocityVerlet,
    cutoff,
    box_size,
)
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
    n = length(simulation_arrays.atom_arrays.r[:])
    # for simpler notation, I shift arrays to new variable names
    r = simulation_arrays.atom_arrays.r[:]
    v = simulation_arrays.atom_arrays.v[:]
    f = simulation_arrays.atom_arrays.f[:]
    m = simulation_arrays.atom_arrays.mass[:]

    dt = integrator.dt
    for i = 1:n
        # update velocities
        v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
        # update positons (apply periodic boundary conditions)
        r[i] = pb!(r[i] .+ v[i] .* dt, box_size[1])
    end
    # update forces
    f = analytical_total_force(simulation_arrays, cutoff, box_size)
    for i = 1:n
        # update velocities
        v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
    end # for loop
end # function

# TODO code in Langevin
struct LangevinIntegrator{T}
    bath_const::T
end

# TODO get langvevin integrator/thermostat working for simulation_arrays
# function apply_integrator!(
#     soa::StructArray,
#     temp,
#     dt,
#     integrator::LangevinIntegrator,
# )
#     b_propagator(soa, dt / 2.0, n)               #! B kick half-step
#     #print("a: ")
#     a_propagator(dt / 2.0, n, soa, boxSize)           #! A drift half-step
#     #print("o: ")
#     o_propagator(soa, dt, n, temp, integrator.bath_const)     #! O random velocities and friction step
#     a_propagator(dt / 2.0, n, soa, boxSize)           #! A drift half-step
#     #print("forces: ")
#     UpdateAllForces!(soa, intraFF, vdwTable, list, point, n, boxSize)
#
#     b_propagator(soa, dt / 2.0, n)
# end


function a_propagator(t::Real, n::Int64, soa, boxSize)
    """Implements the A propagator (drift)"""
    # in: t,soa, boxSize
    # t is typically timestep / 2

    r = Vector{MVector{3,Float64}}(undef, n)
    @inbounds for i = 1:n
        r[i] = soa.r[i] + 0.5 * soa.r[i] * t
    end
    pbc!(soa, r, n, boxSize)
end


function b_propagator(soa, t, n)
    """ Update Velocity """
    @inbounds for i = 1:n
        soa.v[i] = soa.v[i] + soa.f[i] * t
    end
end # b_propagator

function o_propagator(
    soa,
    t::Float64,
    n::Int64,
    temperature::Float64,
    gamma = 2.0,
)
    # in: soa, t, n, gamma, temperature
    # local: c1, c2, c3, c4  # Taylor series coefficients

    c1 = 2.0
    c2 = -2.0
    c3 = 4.0 / 3.0
    c4 = -2.0 / 3.0
    x = gamma * t       # gamma is the friction coefficient
    if (x > 0.0001)  # Use formula
        c = 1 - exp(-2.0 * x)
    else # Use Taylor expansion for low x
        c = x * (c1 + x * (c2 + x * (c3 + x * c4)))
    end
    c = sqrt(c)
    @inbounds for i = 1:n
        soa.v[i] =
            exp(-x) * soa.v[i] + c * rand(Normal(0.0, sqrt(temperature)), 3)
    end

end # o_propagator
