export top_file, specie_location, system_location

top_file = joinpath(
    dirname(pathof(AutoDiff_md)),
    "../",
    "structures",
    "topology_files",
    "N2_lj.top"
    # "Ar.top",
    #"mea_tip3p.top",
) # Ar.top
specie_location = joinpath(
    dirname(pathof(AutoDiff_md)),
    "../",
    "structures",
    "single_molecules",
)  # "mea.pdb",
system_location =
    joinpath(dirname(pathof(AutoDiff_md)), "../", "structures", "whole_systems")
