export Logo, Completion, PrintLine

function PrintLine(s::String, n::Int64) # handy tool for outputting lines
    println(repeat(s, n)) # prints n copies of whatever s is.
end

function Logo()

    PrintLine("=", 70)
    println("                                                                  ")
    println("                   Code Version: Development (0.0)                ")
    println("                      Julia version 1.4.0                         ")
    println("                     Written by Braden Kelly                      ")
    println("                        January, 2021                             ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                       BEGINNING SIMULATION                       ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("        Seriously, you should be holding on to something          ")
    println("                   (Nobody likes a Noob...)                       ")
    println("                                                                  ")
    PrintLine("=", 70)

end

function Completion()

    PrintLine("=", 70)
    println("                                                                  ")
    println("                   Code Version: Development (-1.0)               ")
    println("                      Julia version 1.4.0                         ")
    println("                     Written by Braden Kelly                      ")
    println("                        January, 2021                             ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                       FINISHED SIMULATION                        ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                               .                                  ")
    println("                                                                  ")
    println("                                                                  ")
    println("                                                                  ")
    PrintLine("=", 70)

end
