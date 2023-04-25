using DelimitedFiles
using LightGraphs
using Pkg
using Combinatorics
using DelimitedFiles
using LoopVectorization

function CM()
    degree = reverse(sort([parse(Int, line) for line in eachline("./input/degree.txt")]))
    sitios = 1:length(degree)
    #print(sitios)
	
    lig = falses(length(sitios), length(sitios))
    for i in sitios
        aleatorios = shuffle(sitios[sitios .!= i])
        for j in aleatorios
            if degree[i] == 0
                break
            end
            if degree[j] > 0 && !(@views lig[i, j] || lig[j, i])
                lig[i, j] = true
                degree[i] -= 1
                degree[j] -= 1
                break
            end
        end
    end

    return findall(triu(lig) .== true)
end

function CM100()
    for i in 1:100
        tempo_decorrido = @elapsed CM()
        println(tempo_decorrido)
    end
end

mutable struct Graph
    viz::Matrix{Cint}
    Nodes::Int32
    edges::Int32
end

function SBM(N,seed)
    Random.seed!(seed)
    prob = readdlm("./input/probability.txt",' ')
    P = readdlm("./input/P.txt",' ')

    faixas = round.(Int, P*N)
    if sum(faixas)!= N
        x = rand()*length(faixas)
        x = Int(round(x))
        if x ==0
            x += 1
        end
        faixas[x] += 1
    end

    faixas = shuffle(vcat(repeat([1], faixas[1]), repeat([2], faixas[2]), repeat([3], faixas[3]), repeat([4], faixas[4]), repeat([5], faixas[5])))

    g = SimpleGraph(N)
    for i in 1:N
        sitios = [j for j in i+1:N]
        ligacoes = [prob[faixas[i],faixas[j]] for j in i+1:N]
        r = rand(length(ligacoes))
        ligacoes = [sitios[j] for j in 1:length(ligacoes) if r[j]<= ligacoes[j]]
        for j in ligacoes
            add_edge!(g, i, j)
        end
    end
    file = open("./input/arquivo.txt", "w")
    arquivo = open("./input/faixas.txt", "w")
    for i in 1:N
        vizinhos = neighbors(g,i)
        if(length(vizinhos) > 0)
            if i == 1277
                print(length(vizinhos) > 0)
            end
            vizinhos = [faixas[j] for j in vizinhos]
            cont = [0,0,0,0,0]
            for j in vizinhos
                cont[j] += 1
            end
            cont = join(cont, " ")*"\n"
            write(file, cont)
            s = faixas[i]
            write(arquivo, "$s\n")
        end
    end
    close(file)
    close(arquivo)
end

function Configuration_Model(seed)
    Random.seed!(seed)
    degrees = readdlm("./input/arquivo.txt",' ')
    faixas = readdlm("./input/faixas.txt",' ')

    soma = sum(degrees, dims=2)
    ordem = reverse(sortperm(vec(soma)))
    degrees = degrees[ordem,:]
    faixas = round.(UInt8, faixas)
    faixas = faixas[ordem]

    g = SimpleGraph(length(faixas))
    sitios = 1:length(faixas)
    for i in sitios
        s = shuffle(sitios[sitios .!= i])
        for j in s
            f1 = faixas[i]
            f2 = faixas[j]
            if sum(degrees[i,:]) == 0
                break
            end
            if degrees[i,f2] == 0 || degrees[j,f1] == 0 || has_edge(g, i, j)
                continue
            end
            degrees[i,f2] -= 1
            degrees[j,f1] -= 1
            add_edge!(g, i, j)
        end
    end
    L =maximum([length(neighbors(g,i)) for i in 1:length(faixas)])
    v = Array{Cint}(undef,length(faixas),L+1)
    for i in 1:length(faixas)
        v[i,1] = length(neighbors(g,i))
        v[i,2:length(neighbors(g,i))+1] .= neighbors(g,i)
        v[i,length(neighbors(g,i))+2:end] .= -1
    end
    G = Graph(v,length(faixas),length(collect(edges(g))))
    return G,faixas
end
