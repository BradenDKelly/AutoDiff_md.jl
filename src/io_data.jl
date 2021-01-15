

function ReadCNF(input)
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
