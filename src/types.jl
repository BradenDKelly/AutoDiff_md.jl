abstract type ForceField end
abstract type Thermostat end
abstract type Barostat end
abstract type NeighborList end
abstract type Integrator end
abstract type Gromacs <: ForceField end
abstract type GAFF <: ForceField end
abstract type AngleType <: ForceField end
abstract type DihedralType <: ForceField end
abstract type AtomInfo end
abstract type TypeStructArray end

abstract type molecular_ig end # for virial calculation
abstract type atomic_ig end # for virial calculation
abstract type atomic_ke end # for virial calculation
