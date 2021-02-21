#using StaticArrays
export pb!, apply_molecular_pb!, apply_integrator!, VelocityVerlet

include("types.jl")
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
function pb!(r::SVector{3}, box::SVector{3})

    x = r[1]
    y = r[2]
    z = r[3]

    if x < 0 x += box[1] end
    if x > box[1] x -= box[1] end
    if y < 0 y += box[2] end
    if y > box[2] y -= box[2] end
    if z < 0 z += box[3] end
    if z > box[3] z -= box[3] end

    r = SVector{3}(x, y, z)
    return r
end

""" returns max and min values in an array of SVectors"""
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

function apply_molecular_pb!(sa::SimulationArrays, box_size)
    old_com, new_com = SVector(0, 0, 0), SVector(0, 0, 0)
    sa.molecule_arrays.COM[:] = total_com_update!(sa)
    for i = 1:length(sa.molecule_arrays.COM)
        ia = sa.molecule_arrays[i].firstAtom
        ja = sa.molecule_arrays[i].lastAtom
        #old_com = COM(sa.atom_arrays.r[ia:ja], sa.atom_arrays.mass[ia:ja])
        # TODO figure out this incoherent pointing that Julia is doing.
        old_com = sa.molecule_arrays.COM[i]
        new_com = deepcopy(old_com)
        new_com = pb!(new_com, box_size)
        # recent trust issues have made me put the next line in
        sa.molecule_arrays.COM[i] = new_com
        if new_com != old_com
            for j = ia:ja
                sa.atom_arrays.r[j] = sa.atom_arrays.r[j] + new_com - old_com
            end # for
        end # if
    end # for

    # # for troubleshooting checks that particles are in box
    # minV, maxV = maxmin(sa.molecule_arrays.COM)
    # #println(minV, maxV)
    # if minV < 0.0
    #     println("Shit, particle is less than 0")
    # end
    # if maxV > box_size[1]
    #     println("Shit, particle is outside box")
    # end
    return sa
end

# """Struct with parameters for Velocity Verlet integrator"""
# struct VelocityVerlet{F} <: Integrator
#     dt::F
# end

"""
Velocity-Verlet integrator scheme.

Parameters
----------
simulation_array : SimulationArray
    contains arrays related to atoms, molecules, and their parameters.
cutoff : Float
    interaction cutoff distance
box_size : SVector{3}
    vector of x, y, z box lengths (should all be the same)

Returns
Nothing. Updated in place.
"""
function apply_integrator!(
    sa::SimulationArrays,
    integrator::VelocityVerlet,
    cutoff,
    box_size,
)

    n = length(sa.atom_arrays.r[:])
    # for simpler notation, I shift arrays to new variable names
    r = sa.atom_arrays.r[:]
    v = sa.atom_arrays.v[:]
    f = sa.atom_arrays.f[:]
    m = sa.atom_arrays.mass[:]

    dt = integrator.dt
    for i = 1:n
        # update velocities
        v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
        # update positons (apply periodic boundary conditions)
        #r[i] = pb!(r[i] .+ v[i] .* dt, box_size)
        r[i] = r[i] .+ v[i] .* dt
    end
    sa.atom_arrays.r[:] = r
    # TODO remove mirror image in intramolecular calculations
    apply_molecular_pb!(sa, box_size)
    # TODO figure out why v doesn't point to sa.atom_arrays.v and vice versa
    #sa.atom_arrays.v[:] = v
    # update forces
    f = analytical_total_force(
        sa,
        cutoff,
        box_size,
        sa.neighborlist.point[:],
        sa.neighborlist.list[:],
    )
    for i = 1:n
        # update velocities
        v[i] = v[i] .+ 0.5 .* dt .* f[i] ./ m[i]
    end # for loop
    # TODO figure out why v doesn't point to sa.atom_arrays.v and vice versa
    sa.atom_arrays.v[:] = v
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

"""Implements the A propagator (drift)"""
function a_propagator(t::Real, n::Int64, soa, boxSize)
    # in: t,soa, boxSize
    # t is typically timestep / 2

    r = Vector{MVector{3,Float64}}(undef, n)
    @inbounds for i = 1:n
        r[i] = soa.r[i] + 0.5 * soa.r[i] * t
    end
    pbc!(soa, r, n, boxSize)
end

""" Update Velocity """
function b_propagator(soa, t, n)
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
