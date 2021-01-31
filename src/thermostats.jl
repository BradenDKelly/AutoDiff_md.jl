export
    AndersenThermostat,
    apply_thermostat

struct AndersenThermostat{T} <: Thermostat
    bath_couple::T
end

function apply_thermostat!(v, m, temp, dt, thermostat::AndersenThermostat)
    n = length(v)
    for i=1:n
        # update velocities
        if rand() < thermostat.bath_couple * dt
            # Andersen Thermostat
            v[i] = velocity(m[i], temp)
        end # if statement
    end # for loop
end

struct NoThermostat  <: Thermostat end

function apply_thermostat!(v, m, temp, dt, thermostat::NoThermostat)
    return v
end
