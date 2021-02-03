export Read_gromacs, ReadCNF, PrintPDB_argon

"""
Reads initial structures and generates position and velocity arrays

Parameters
-----------
input : string
    name of file to open

Notes
----------
some files store position then quaternion elements (4 elements to a quat)
in this case we just assume 3 of them as our initial velocity even though
this clearly makes no sense, we just need a start. We could override
with Maxwell-Boltzmann if we wanted.

hence, "e" is actually for quaternion.
"""
function ReadCNF(input)
    r = []
    e = []
    i = 0
    box1 = 0.0
    #box=9.42953251
    open(input) do file
        for line in eachline(file)
            i += 1

            if i == 2 #length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = parse(Float64, strip(line)) # 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
                #println("hardcoded box at line 257 in ReadCNF: ", box1)
            end

            if i >= 3 #length(split(line)) > 4

                lin = split(line)

                push!(
                    r,
                    [
                        parse(Float64, strip(lin[1])),
                        parse(Float64, strip(lin[2])),
                        parse(Float64, strip(lin[3])),
                    ],
                )

                push!(
                    e,
                    [
                        parse(Float64, strip(lin[4])),
                        parse(Float64, strip(lin[5])),
                        parse(Float64, strip(lin[6])),
                        #parse(Float64, strip(lin[7])), # commented out Jan 14 2021
                    ],
                )
            end
        end
    end
    r = [SVector{3,Float64}(r[i]...) for i = 1:length(r)]
    return r, e, box1
end

"""
Reads initial coordinates and velocity from a gromacs file

Parameters
-----------
input : string
    name of file to open

Notes
----------
some files store position then quaternion elements (4 elements to a quat)
in this case we just assume 3 of them as our initial velocity even though
this clearly makes no sense, we just need a start. We could override
with Maxwell-Boltzmann if we wanted.

hence, "e" is actually for quaternion.
"""
function Read_gromacs(input)

    r = []
    v = []
    i = 0
    box1 = 0.0
    #box=9.42953251
    position = true
    velocity = false
    box = false
    open(input) do file
        for line in eachline(file)
            if i > 5 && occursin("end", lowercase(line))
                position = false
            elseif i > 5 && occursin("velocity", lowercase(line))
                velocity = true
            end
            i += 1

            if i >= 5 && position#length(split(line)) > 4

                lin = split(line)

                push!(
                    r,
                    [
                        parse(Float64, strip(lin[5])),
                        parse(Float64, strip(lin[6])),
                        parse(Float64, strip(lin[7])),
                    ],
                )
            elseif i > 5 && velocity && length(split(line)) > 3
                lin = split(line)

                push!(
                    v,
                    [
                        parse(Float64, strip(lin[5])),
                        parse(Float64, strip(lin[6])),
                        parse(Float64, strip(lin[7])),
                    ],
                )
            end # if
            if occursin("box", lowercase(line))
                box = true
            end
            if i > 5 && box && length(split(line)) == 3
                lin = split(line)
                # println(line)
                box1 = parse(Float64, strip(lin[1]))
                box = false
            end
        end # for
    end # open
    r = [SVector{3,Float64}(r[i]...) for i = 1:length(r)]
    v = [SVector{3,Float64}(v[i]...) for i = 1:length(v)]
    #@assert length(r) == length(v)
    return r, v, box1
end

""" Print a pdb file for this snapshot. This is for monatomic system"""
function PrintPDB_argon(r::Vector, boxSize, step = 1, filename = "pdbOutput")

    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf(
            "%-7s %7.3f %7.3f %7.3f %30s \n",
            "CRYST1",
            10.0 * boxSize[1],
            10.0 * boxSize[2],
            10.0 * boxSize[3],
            "90.00  90.00  90.00 P 1           1"
        )
        write(file, line)

        for (i, coord) in enumerate(r)
            atomName = "Ar" #systemTop.molParams[soa.molType[i]].atoms[soa.atype[i]].atomnm
            molName = "Ar" #systemTop.molParams[soa.molType[i]].atoms[soa.atype[i]].resnm

            line = @sprintf(
                "%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
                "ATOM",
                i,
                atomName,
                molName,
                i,
                " ",
                10.0 * coord[1],
                10.0 * coord[2],
                10.0 * coord[3],
                1.00,
                0.00
            )
            write(file, line)
        end
    end
end
