using Random
using LinearAlgebra
function CM()
    degree = reverse(sort([parse(Int, line) for line in eachline("./degree.txt")]))
    sitios = 1:length(degree)
    #print(sitios)
	
    lig = falses(length(sitios), length(sitios))

    for i in sitios
        println(i)
        while degree[i] > 0
            println(i," ",degree[i])
            aleatorios = shuffle(sitios[sitios .!= i])
            for j in aleatorios
                if degree[j] > 0 && !(@views lig[i, j] || lig[j, i])
                    lig[i, j] = true
                    degree[i] -= 1
                    degree[j] -= 1
                    break
                end
            end
            degree[i] == 0 && break
        end
    end

    return findall(triu(lig) .== true)
end

tempo_decorrido = @elapsed CM()
println("Tempo decorrido: ", tempo_decorrido, " segundos")
