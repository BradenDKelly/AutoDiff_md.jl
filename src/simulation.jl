using StaticArrays

include("forces.jl")
include("energies.jl")
include("integrators.jl")
include("thermo_functions.jl")


function simulate(r, v, box_size, temp, dt, nsteps)
    """
    The main loop for generating the molecular dynamics trajectory

    Parameters
    ----------
    r : Vector{SVector{3}}
        vector of all atom coordinates
    box_size : SVector{3}
        Vector of x, y, z box lengths
    temp : Float64
        set temperature (not used, we don't use Maxwell-Boltzmann for initial
        velocities, we read them from file)
    dt : Float64
        time step (approximately 0.005 for reduced units)
    nsteps : Int64
        number of steps to take
    v : Vector{SVectors}
        Starting velocites for each atom - read from some file

    Returns
    ----------
    Currently just prints samples, but if it was actually to be used, it would
    return sample results
    """
    n = length(r)
    v = initial_velocity(v)
    println("Velocity center of mass: ", sum(v))
    println("Initial temperature is $(dot(v,v)/(3 * n - 3)) and A&T is 0.223094")
    temp_stat=[]
    press_stat=[]
    ener_stat = []
    t = 0
    count = 0
    f = analytical_total_force(r, box_size)
    vol = box_size[1]^3
    println("Initial pressure: ", n/vol + virial(r, box_size) / vol / 3, " A&T: ", 0.32 + (-603)/box_size[1]^3)


    for i=1:nsteps
        # Calculate forces on all atoms
        f = analytical_total_force(r, box_size)
        # Solve Newtons equations of motion (we use Velocity Verlet)
        integrator!(r, v, f, dt, box_size)
        count += 1
        # Sample some things
        if count == 1000
            println("Current step is $i")
            push!(temp_stat, dot(v,v)/(3 * n - 3))
            KE = kinetic_energy(v)
            tmp = 2 * KE / (3 * n - 3)
            tot_vir = virial(r, box_size)
            press_f = press_full(tot_vir, n, vol, tmp)
            push!(press_stat, press_f)
            push!(ener_stat, KE + total_energy(r, box_size)[1])
            count = 0
        end
        # step time forward
        t += dt
    end
    println("average temperature was: ", mean(temp_stat))
    println("average pressure was: ", mean(press_stat), " std is: ", std(press_stat))
    println("average energy was: ", mean(ener_stat))
    println("average energy pp was: ", mean(ener_stat)/n, " std is: ", std(ener_stat)/n)
end
