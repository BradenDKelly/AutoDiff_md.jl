#using StaticArrays
export
    simulate

include("forces.jl")
include("energies.jl")
include("integrators.jl")
include("thermo_functions.jl")


function simulate(r, v, m, eps, sig, box_size, temp, dt, nsteps, cutoff)
    """
    The main loop for generating the molecular dynamics trajectory

    Parameters
    ----------
    r : Vector{SVector{3}}
        vector of all atom coordinates
    v : Vector{SVectors}
        Starting velocites for each atom - read from some file
    m : Vector
        masses of all atoms
    eps : Vector
        epsilon parameter of all atoms
    sig : Vector
        sigma parameters of all atoms
    box_size : SVector{3}
        Vector of x, y, z box lengths
    temp : Float64
        set temperature (not used, we don't use Maxwell-Boltzmann for initial
        velocities, we read them from file)
    dt : Float64
        time step (approximately 0.005 for reduced units)
    nsteps : Int64
        number of steps to take
    cutoff : Float
        interaction cutoff distance


    Returns
    ----------
    Currently just prints samples, but if it was actually to be used, it would
    return sample results
    """
    n = length(r)
    v = initial_velocity(m, v)
    println("Momentum center of mass: ", sum(m .* v))
    #println("Initial temperature is $(dot(v,v)/(3 * n - 3)) and A&T is 0.223094")
    temp_stat=[]
    press_stat=[]
    ener_stat = []
    t = 0
    f = analytical_total_force(r, eps, sig, cutoff, box_size)
    vol = box_size[1]^3
    println("Initial pressure: ", n/vol + virial(r, eps, sig, cutoff, box_size) / vol / 3, " A&T: ", 0.32 + (-603)/box_size[1]^3)
    println("Initial KE $(kinetic_energy(v, m))")
    println("Initial temperature: $(2 * kinetic_energy(v, m) / (3 * n - 3))")

    for i=1:nsteps
        # Calculate forces on all atoms
        f = analytical_total_force(r, eps, sig, cutoff, box_size)
        # Solve Newtons equations of motion (we use Velocity Verlet)
        integrator!(r, v, f, m, dt, eps, sig, cutoff, box_size)
        # Sample some things
        if i % 100 == 0
            KE = kinetic_energy(v, m)
            tmp = 2 * KE / (3 * n - 3)
            push!(temp_stat, tmp)
            tot_vir = virial(r, eps, sig, cutoff, box_size)
            press_f = press_full(tot_vir, n, vol, tmp)
            push!(press_stat, press_f)
            push!(ener_stat, KE + total_energy(r, eps, sig, cutoff, box_size)[1])
            if i % 1000 == 0
                println("Current step is $i, temperature is $tmp, pressure is $press_f total energy is $(KE + total_energy(r, eps, sig, cutoff, box_size)[1])
                ")
            end
        end


        # step time forward
        t += dt
    end
    println("average temperature was: ", mean(temp_stat))
    println("average pressure was: ", mean(press_stat), " std is: ", std(press_stat))
    println("average energy was: ", mean(ener_stat))
    println("average energy per particle (E/n) was: ", mean(ener_stat)/n, " std is: ", std(ener_stat)/n)
    display(plot(1:length(temp_stat), temp_stat))
    #plim = box_size[1]
    #display(plot3d(r, xlim = (0, plim), ylim = (0, plim), zlim = (0, plim),title = "Simulation",marker = 2))
end
