using Random
using LinearAlgebra
using Pkg
using Libdl
using LightGraphs
using Statistics
using PackageCompiler
include("rede.jl")
include("infect.jl")
using CxxWrap
# criar um grafo aleatório usando a distribuição de grau dada

#tempo_decorrido = @elapsed CM()
#tempo_decorrido = @elapsed roda()
#println("Tempo decorrido:",tempo_decorrido/100) =#
#SBM(2029,45655)


#tempo_decorrido = @elapsed generate_infect(1)
#println("Tempo decorrido:",tempo_decorrido)




# Caminho para o arquivo funcoes.jl
# Cria a imagem do sistema com as funções do arquivo funcoes.jl

# Set the path to funcoes.jl
#path_to_mylib = joinpath(dirname(pathof(mylib)), "..")
function minha_funcao(viz::Vector{Ptr{Cint}},y::Cint)
    g = SimpleGraph(y)
    c = 0
    println(unsafe_load(viz)[1,1]," ", y)
end
ptr = @cfunction(minha_funcao, Cvoid, (Vector{Ptr{Cint}},Cint))
lib = dlopen("./main.so")
func = dlsym(lib,"main")
ccall(func,Cvoid, (Ptr{Nothing},),ptr)