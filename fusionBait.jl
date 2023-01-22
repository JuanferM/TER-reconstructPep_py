function readFile(path, baits)
    open(path) do file
        k::Int64 = parse(Int64, readline(file))
        for i in 1:k
            bait::String, n::Int64 = split(readline(file))
            baits[bait] = []
            for j in 1:n
                push!(baits[bait], readline(file))
            end
        end
    end
end

function main()
    for baitCollection in ARGS
        baits = Dict{String, Vector{String}}()
        readFile(baitCollection, baits)

    end
end

main()
