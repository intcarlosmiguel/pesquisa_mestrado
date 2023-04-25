using Random
using LinearAlgebra
using Pkg
using Libdl
using LightGraphs
using Statistics
include("rede.jl")
include("infect.jl")

#= function CM100()
    for i in 1:100
        tempo_decorrido = @elapsed CM()
        println(tempo_decorrido)
    end
end

function Configuration_Model(degrees)
    for i in 100:150
        println(i)
        g = LightGraphs.SimpleGraphs.random_configuration_model(length(degrees),degrees,seed = i)
        graus = degree(g)
    end

end =#


# criar um grafo aleatório usando a distribuição de grau dada

#tempo_decorrido = @elapsed CM()
#tempo_decorrido = @elapsed roda()
#println("Tempo decorrido:",tempo_decorrido/100) =#
#SBM(2029,45655)

G,faixas =  Configuration_Model(455)
const libm = Libdl.find_library("./main.so")
faixas = convert(Vector{Cint},faixas)
faixas = faixas .-1
ccall((:generate_infect, "./main.so"), Cvoid, (Graph,Ptr{Cint}, Cint,Cint,Cfloat), G,faixas, 42,1,0.0)