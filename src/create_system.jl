using StaticArrays

function initial_velocity(e::Vector)
    """ This assumes we are being given a set of velocities.
    We shift the velocities so that the center of mass velocity is 0

    Parameters
    ----------
    e : Vector
        vector of SVectors, the list of each atoms starting velocities as read
        from some file

    Returns
    ----------
    v : Vector{SVector{3}}
        vector of velocites that have center of mass velocity of zero
    """
    n = length(e)
    v = [SVector{3,Float64}(e[i]...) for i = 1:length(e)]
    sumv = sum(v)./n
    return [(v[i].-sumv) for i=1:n]
end

function InitCubicGrid(n::Int, rho::Real)

    #!------------------------------------------------------------------------
    #! Created by Braden Kelly
    #!------------------------------------------------------------------------
    #! Creates an initial configuration
    #!------------------------------------------------------------------------
    #! input:  n       number of particles
    #!         rho     density
    #! output: coords  coordinates of n particles
    #!         L       box length
    #!------------------------------------------------------------------------

    #! Calculate box length (L)
    L = (n / rho)^(1.0 / 3.0)

    #! Calculate the lowest perfect cube that will contain all of the particles
    nCube = 2

    while (nCube^3 < n)
        nCube = nCube + 1
    end
    coords = Array{Float64,2}(undef, 3, n)
    # initial position of particle 1
    posit = [0, 0, 0]

    #!begin assigning particle positions
    for i = 1:n
        coords[:, i] = (posit + [0.01, 0.01, 0.01]) * (L / nCube)
        #! Advancing the index (posit)
        posit[1] = posit[1] + 1
        if (posit[1] == nCube)
            posit[1] = 0
            posit[2] = posit[2] + 1
            if (posit[2] == nCube)
                posit[2] = 0
                posit[3] = posit[3] + 1
            end
        end
    end
    return [
        SVector{3}(coords[1, i], coords[2, i], coords[3, i])
        for i = 1:size(coords, 2)
    ], L
end #function
