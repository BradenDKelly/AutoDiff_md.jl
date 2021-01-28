export
    andersen!

function andersen!(v, m, temp, dt, thermo_couple)
    n = length(v)
    for i=1:n
        # update velocities
        if rand() < thermo_couple * dt
            # Andersen Thermostat
            v[i] = velocity(m[i], temp)
        end # if statement
    end # for loop
end
