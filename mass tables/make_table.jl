using Printf
using Base.Threads

code = ["", "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L",
        "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U"]

mono = Dict{String, Float64}("A" => 71.037113805,
        "R" => 156.10111105,
        "N" => 114.04292747,
        "D" => 115.026943065,
        "C" => 103.009184505,
        "E" => 129.042593135,
        "Q" => 128.05857754,
        "G" => 57.021463735,
        "H" => 137.058911875,
        "I" => 113.084064015,
        "L" => 113.084064015,
        "K" => 128.09496305,
        "M" => 131.040484645,
        "F" => 147.068413945,
        "P" => 97.052763875,
        "S" => 87.032028435,
        "T" => 101.047678505,
        "W" => 186.07931298,
        "Y" => 163.063328575,
        "V" => 99.068413945,
        "U" => 150.953633405)

comb = Dict()
mass = Dict()

println("Computing...")
counter = 0
a = time()
for i in code
    for j in code
        for k in code
            for l in code
                for m in code
                    for n in code
                        for o in code
                            @threads for p in code
                                global counter += 1
                                word = join(sort(collect(i*j*k*l*m*n*o*p)))
                                if word != "" && !haskey(comb,word)
                                    comb[word] = true
                                    s = 0
                                    for c in word
                                        if c != ""
                                            s += mono[string(c)]
                                        end
                                    end
                                    s = trunc(s, digits=2)
                                    if haskey(mass,s)
                                        push!(mass[s],word)
                                    else
                                        mass[s] = [word]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

println("Done.\nWriting to csv...")
open("table.csv", "w") do f
    write(f, string(length(mass))*"\n")
    # sort!(collect(mass), by = x->x[1])
    for (k, v) in mass
        write(f, string(k)*",")
        for i in v
            write(f, i*",")
        end
        write(f, "\n")
    end
end
println("Done.")
b = time()
@printf("elapsed time is %.2f sec.", b-a)
