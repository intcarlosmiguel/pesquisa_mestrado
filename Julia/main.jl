using Random
using LinearAlgebra
using Pkg
using Libdl
using LightGraphs
using Statistics
include("rede.jl")
include("infect.jl")


# criar um grafo aleatório usando a distribuição de grau dada

#tempo_decorrido = @elapsed CM()
#tempo_decorrido = @elapsed roda()
#println("Tempo decorrido:",tempo_decorrido/100) =#
#SBM(2029,45655)


tempo_decorrido = @elapsed generate_infect(1)
println("Tempo decorrido:",tempo_decorrido)