#using StaticArrays
export simulate!

include("forces.jl")
include("energies.jl")
include("integrators.jl")
include("thermo_functions.jl")

function simulate!(
    simulation_arrays::SimulationArrays,
    simulation_controls::SimulationControls,
    box_size::SVector,
    nsteps,
    cutoff,
    props
)
    """
    The main loop for generating the molecular dynamics trajectory

    Parameters
    ----------
    simulation_array : SimulationArray
        atom_arrays::StructArray
            molNum::I
            molType::I
            atype::I
            mass::F
            r::SVector{3,F}
            v::SVector{3,F}
            f::SVector{3,F}
            qq::Float64
        molecule_arrays::StructArray
        intraFF::IntraForceField
            bonds::Vector{Bonds}
            angles::Vector{Angles}
            dihedrals::Vector{Dihedrals}
        vdwTable::Tables
            eps : table of lj epsilon parameters
            sig : table of lj sigma parameters
        nonbonded_matrix::Array
        scaled_pairs::Vector
        numbers::Numbers
            atoms::I
            molecules::I
            atomTypes::I
            molTypes::I
            charges::I
        neighborlist::NeighborList
    simulation_controls : SimulationControls
        integrator::VelocityVerlet
        thermostat::Thermostat
        barostat::Barostat
    box_size : SVector{3}
        Vector of x, y, z box lengths
    nsteps : Int64
        number of steps to take
    cutoff : Float
        interaction cutoff distance
    props : Struct
        struct with properties such as temp, press, energy etc

    Returns
    ----------
    Currently just prints samples, but if it was actually to be used, it would
    return sample results
    """
    n = length(simulation_arrays.atom_arrays.r)
    simulation_arrays.atom_arrays.v[:] =
        [velocity(
            simulation_arrays.atom_arrays.mass[i],
            simulation_controls.thermostat.set_temp
        ) for i=1:n]

    simulation_arrays.atom_arrays.v[:] = initial_velocity(
        simulation_arrays.atom_arrays.mass,
        simulation_arrays.atom_arrays.v
    )
    vol = volume(box_size)
    temp_stat = []
    press_stat = []
    ener_stat = []
    KE_stat = []
    PE_stat = []
    t = 0

    for i = 1:nsteps

        # Calculate forces on all atoms
        f = analytical_total_force(simulation_arrays, cutoff, box_size)
        #f = analytical_total_force(r, atype, vdwTable, cutoff, box_size)
        simulation_arrays.atom_arrays.f[:] = f
        # Solve Newtons equations of motion (use vv and Andersen thermostat)
        apply_integrator!(
            simulation_arrays,
            simulation_controls.integrator,
            cutoff,
            box_size
        )
        # thermostat
        apply_thermostat!(
            simulation_arrays.atom_arrays.v,
            simulation_arrays.atom_arrays.mass,
            simulation_controls.integrator.dt,
            simulation_controls.thermostat
        )

        # barostat
        if simulation_controls.barostat.use

            if i % simulation_controls.barostat.pcouple == 0

                r, box_size, cutoff = apply_barostat!(
                    simulation_arrays,
                    simulation_controls.barostat,
                    box_size,
                    cutoff,
                    simulation_controls.thermostat.set_temp
                )
                simulation_arrays.atom_arrays.r[:] = r
                vol = volume(box_size)
            end
        end
        # Sample some things
        if i % 100 == 0
            v = simulation_arrays.atom_arrays.v[:]
            m = simulation_arrays.atom_arrays.mass[:]
            r = simulation_arrays.atom_arrays.r[:]
            KE = kinetic_energy(v, m)
            push!(KE_stat, KE)
            tmp = 2 * KE / (3 * n - 3) / 0.0083144621
            push!(temp_stat, tmp)
            tot_vir = virial(simulation_arrays, cutoff, box_size)
            press_f = press_full(tot_vir, n, vol, tmp)
            #PE = total_energy(r, 0.9906, 0.3405, cutoff, box_size)[1]
            push!(press_stat, press_f)
            #push!(ener_stat, KE + PE)
            #push!(PE_stat, PE)
            if i % 1000 == 0
                @printf(
                    "Current step is %d, temp is %.3f, KE is %.3f, press is %.3f, box is %.3f\n",
                    i,
                    tmp,
                    KE,
                    #PE,
                    press_f,
                    #KE + PE,
                    box_size[1]
                )
                #println("coords $(r[1]), $(r[250]), $(r[500]), $(r[750]), $(r[1000])")
            end
        end

        # if i % 1000 == 0
        #     PrintPDB_argon(r, box_size, i)
        # end


        # step time forward
        t += simulation_controls.integrator.dt
    end
    println("average temperature was: ", mean(temp_stat))
    println(
        "average pressure was: ",
        mean(press_stat),
        " std is: ",
        std(press_stat),
    )
    # println("average energy was: ", mean(ener_stat))
    println("average KE was: ", mean(KE_stat))
    println("average PE was: ", mean(PE_stat))
    println(
        "average energy per particle (E/n) was: ",
        mean(ener_stat) / n,
        " std is: ",
        std(ener_stat) / n,
    )
    display(plot(
        1:length(temp_stat),
        temp_stat,
        title = "Temperature vs steps(100)",
        label = "Temperature",
        xlabel = "Iteration (100's)",
        ylabel = "Kinetic Temperature, K",
    ))
    display(plot(
        1:length(KE_stat),
        KE_stat,
        title = "KE vs steps(100)",
        label = "KE",
        xlabel = "Iteration (100's)",
        ylabel = "Kinetic Energy kJ/mol",
    ))
    display(plot(
        1:length(PE_stat),
        PE_stat,
        title = "PE vs steps(100)",
        label = "PE",
        xlabel = "Iteration (100's)",
        ylabel = "Potential Energy, kJ/mol",
    ))
    props.kinetic_temp = mean(temp_stat)
    props.pressure = mean(press_stat)
    props.total_energy = mean(ener_stat)
    #plim = box_size[1]
    #display(plot3d(r, xlim = (0, plim), ylim = (0, plim), zlim = (0, plim),title = "Simulation",marker = 2))
end

function simulate!(
    r,
    v,
    m,
    eps,
    sig,
    box_size,
    temp,
    press,
    dt,
    nsteps,
    cutoff,
    props,
    pcouple,
)
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
    press : Float64
        set (desired) pressure for the npt ensemble
    dt : Float64
        time step (approximately 0.005 for reduced units)
    nsteps : Int64
        number of steps to take
    cutoff : Float
        interaction cutoff distance
    props : Struct
        struct with properties such as temp, press, energy etc


    Returns
    ----------
    Currently just prints samples, but if it was actually to be used, it would
    return sample results
    """
    n = length(r)
    v = initial_velocity(m, v)
    println("Momentum center of mass: ", sum(m .* v))
    #println("Initial temperature is $(dot(v,v)/(3 * n - 3)) and A&T is 0.223094")
    temp_stat = []
    press_stat = []
    ener_stat = []
    KE_stat = []
    PE_stat = []
    t = 0
    f = analytical_total_force(r, eps, sig, cutoff, box_size)
    vol = box_size[1]^3
    println(
        "Initial pressure: ",
        n / vol + virial(r, eps, sig, cutoff, box_size) / vol / 3,
        " A&T: ",
        0.32 + (-603) / box_size[1]^3,
    )
    println("Initial KE $(kinetic_energy(v, m))")
    println("Initial temperature: $(2 * kinetic_energy(v, m) / (3 * n - 3) / 0.0083144621)")
    vmax = log(box_size[1]^3 * 0.01)
    println("vmax $vmax exp(vmax) $(exp(vmax))")
    thermo = true
    thermostat = AndersenThermostat(1.0)
    barostat = true
    v_attempt = 0
    v_accept = 0
    for i = 1:nsteps
        # Calculate forces on all atoms
        f = analytical_total_force(r, eps, sig, cutoff, box_size)

        # Solve Newtons equations of motion (we use Velocity Verlet)
        #integrator!(r, v, f, m, dt, eps, sig, cutoff, box_size)

        # Solve Newtons equations of motion (use vv and Andersen thermostat)
        integrator!(r, v, f, m, dt, eps, sig, cutoff, box_size, temp)
        # thermostat
        if thermo
            apply_thermostat!(v, m, temp, dt, thermostat)
        end

        # barostat
        if barostat && i % pcouple == 0
            v_attempt += 1
            box_old = deepcopy(box_size[1])
            r_old = deepcopy(r[1])
            r, box_size, cutoff = volume_change_lj_MC!(
                r,
                box_size,
                cutoff,
                eps,
                sig,
                temp,
                press,
                vmax,
            )
            if box_size[1] != box_old
                v_accept += 1
            end
            vol = box_size[1] * box_size[2] * box_size[3]
        end
        # Sample some things
        if i % 100 == 0
            KE = kinetic_energy(v, m)
            push!(KE_stat, KE)
            tmp = 2 * KE / (3 * n - 3) / 0.0083144621
            push!(temp_stat, tmp)
            tot_vir = virial(r, eps, sig, cutoff, box_size)
            press_f = press_full(tot_vir, n, vol, tmp)
            PE = total_energy(r, eps, sig, cutoff, box_size)[1]
            push!(press_stat, press_f)
            push!(ener_stat, KE + PE)
            push!(PE_stat, PE)
            if i % 1000 == 0
                @printf(
                    "Current step is %d, temp is %.3f, KE is %.3f, PE is %.3f, press is %.3f, total energy is %.3f box is %.3f acceptance ratio %.3f\n",
                    i,
                    tmp,
                    KE,
                    PE,
                    press_f,
                    KE + PE,
                    box_size[1],
                    v_accept / v_attempt
                )
                #println("coords $(r[1]), $(r[250]), $(r[500]), $(r[750]), $(r[1000])")
            end
        end

        if i % 1000 == 0
            PrintPDB_argon(r, box_size, i)
        end


        # step time forward
        t += dt
    end
    println("average temperature was: ", mean(temp_stat))
    println(
        "average pressure was: ",
        mean(press_stat),
        " std is: ",
        std(press_stat),
    )
    println("average energy was: ", mean(ener_stat))
    println("average KE was: ", mean(KE_stat))
    println("average PE was: ", mean(PE_stat))
    println(
        "average energy per particle (E/n) was: ",
        mean(ener_stat) / n,
        " std is: ",
        std(ener_stat) / n,
    )
    display(plot(
        1:length(temp_stat),
        temp_stat,
        title = "Temperature vs steps(100)",
        label = "Temperature",
        xlabel = "Iteration (100's)",
        ylabel = "Kinetic Temperature, K",
    ))
    display(plot(
        1:length(KE_stat),
        KE_stat,
        title = "KE vs steps(100)",
        label = "KE",
        xlabel = "Iteration (100's)",
        ylabel = "Kinetic Energy kJ/mol",
    ))
    display(plot(
        1:length(PE_stat),
        PE_stat,
        title = "PE vs steps(100)",
        label = "PE",
        xlabel = "Iteration (100's)",
        ylabel = "Potential Energy, kJ/mol",
    ))
    props.kinetic_temp = mean(temp_stat)
    props.pressure = mean(press_stat)
    props.total_energy = mean(ener_stat)
    #plim = box_size[1]
    #display(plot3d(r, xlim = (0, plim), ylim = (0, plim), zlim = (0, plim),title = "Simulation",marker = 2))
end
