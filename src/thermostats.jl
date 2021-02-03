export AndersenThermostat, apply_thermostat

# TODO add more thermostats i.e., Berendsen and Nose-Hoover chains
"""Struct with parameters for Andersen thermostat"""
struct AndersenThermostat{T} <: Thermostat
    bath_couple::T
    set_temp::T
end

"""Applies the passed thermostat"""
function apply_thermostat!(v, m, dt, thermostat::AndersenThermostat)
    n = length(v)
    for i = 1:n
        # update velocities
        if rand() < thermostat.bath_couple * dt
            # Andersen Thermostat
            v[i] = velocity(m[i], thermostat.set_temp)
        end # if statement
    end # for loop
end

"""Struct for handling no thermostat"""
struct NoThermostat <: Thermostat end

"""Applies no thermostat"""
function apply_thermostat!(v, m, temp, dt, thermostat::NoThermostat)
    return v
end
