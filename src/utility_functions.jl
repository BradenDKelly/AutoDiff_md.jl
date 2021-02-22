export FindMolType,
    FindAtomInMol,
    FindNumericAtomType,
    COM,
    Shift_COM_to_Zero!,
    total_com_update!

include("structs.jl")

# TODO this is not robust, if different molecules have the same atomtypes in them...
function FindMolType(type::String, moleculeList::Array)

    for (i, molecule) in enumerate(moleculeList)
        for (j, atom) in enumerate(molecule.resnm)
            if lowercase(atom) == lowercase(type)
                return i
            end # if
        end # for
    end # for
end # function

function FindAtomInMol(atomType::String, molType::Int64, systemTop)
    for (i, item) in enumerate(systemTop.molParams[molType].atoms)
        if lowercase(item.atomnm) == lowercase(atomType)
            return i
        end
    end
end

function FindNumericAtomType(molType::Int64, atomType::Int64, systemTop)
    name = systemTop.molParams[molType].atoms[atomType].type
    for (i, item) in enumerate(systemTop.atomTypes)
        if lowercase(item.name) == lowercase(name)
            return i
        end
    end
end

function COM(atoms::Vector, masses)
    totalMass = sum(masses)
    numerator = sum(atoms .* masses)
    return numerator ./ totalMass
end

#"Calculates the Center of Mass of a molecule"
function Center_of_Mass(atom_coords, mol_mass) #result (COM)
    #!============================================================================================================

    denominator = 0.0       #! sum of atom masses
    numerator = zeros(3)         #! weighted sum of atom masses and positions

    for i = 1:length(mol_mass)
        numerator .+= atom_coords[i] * mol_mass[i]
        denominator += mol_mass[i]

    end

    return SVector(numerator / denominator...)

end #Center_of_Mass
#!=============================

#!============================================================================================================
function Shift_COM_to_Zero!(atm_coords, COM)

    #!============================================================================================================

    # real, dimension(3), intent(in)                         			:: COM
    # integer, intent(in)                                    			:: Atoms_in_molecule
    # real, intent(inout), dimension(3,atoms_in_molecule)    			:: atm_coords							!atm = atom, qq = charge
    # integer                                                			:: i
    #!-----------------------------------------------------------------------------------------------------------

    #! This is meant to take the body fixed coordinates of the atoms for a molecule and shift them so that the center of mass is at zero. This is necessary for proper rotation of the molecule.
    #! This should only be called at the start of the run after getting *.pdb file coordinates.

    for i = 1:length(atm_coords)
        atm_coords[i] = SVector(atm_coords[i] .- COM[:])

    end

    return nothing

end #Shift_COM_to_Zero

#"""Update all centers-of-mass"""
function total_com_update!(sa::SimulationArrays)
    for (indx, com) in enumerate(sa.molecule_arrays.COM)
        i = sa.molecule_arrays[indx].firstAtom
        j = sa.molecule_arrays[indx].lastAtom
        #println(i, j, typeof(i), typeof(j))
        #println(Center_of_Mass(sa.atom_arrays.r[i:j], sa.atom_arrays.mass[i:j]))
        sa.molecule_arrays.COM[indx] = COM(sa.atom_arrays.r[i:j], sa.atom_arrays.mass[i:j])
    end
    return sa.molecule_arrays.COM
end
