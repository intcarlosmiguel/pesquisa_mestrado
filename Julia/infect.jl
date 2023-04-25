include("rede.jl")
using Pkg
using StaticArrays: @SVector

beta1::Float64 = 0.5;
beta2::Float64 = 0.41;
sigma::Float64 = 1/5.1;
gamma_A::Float64 = 1/7;
phi::Float64 = 0.025;
gamma_I::Float64 = 1/7;
gamma_H::Float64 = 1/14;
delta::Float64 = 1/14;
recupera::Float64 = 1/40;

function calc_estagio(site,g,estagio,vacinado,faixas,prob_estagio,morte,hospitalizacao)
    vacinado_site::Float64 = vacinado[site] == 1 ? 0.058 : 1


    if estagio[site] == 0
        beta  = mapreduce(i -> (estagio[i] == 1 ? beta1*(vacinado[i] == 1 ?  0.058 : 1) : estagio[i] == 2 ? beta2*(vacinado[i] == 1 ?  0.058 : 1)  : 0),+, neighbors(g, site))
        return beta*vacinado_site
    elseif estagio[site] == 1

        return sigma

    elseif estagio[site] == 2

        return gamma_A

    elseif estagio[site] == 3

        hs = hospitalizacao[faixas[site]]
        if vacinado[site] == 1
             hs *= 0.034
        end
        if(prob_estagio[site]<hs) 
            return phi
        else 
            return gamma_I
        end

    elseif estagio[site] == 4

        m = morte[faixas[site]]
        if vacinado[site] == 1
             m *= 0.034
        end
        if(prob_estagio[site]<m) 
            return delta
        else 
            return gamma_H
        end

    elseif estagio[site] == 5
        return recupera
    else
        return 0
    end

end

function exponentialRand(lambda)
    return lambda != 0 ? -log(1 - rand()) / lambda : 0
end


function infect(seed,S0,E0,infect_time,quant)
    Random.seed!(seed)
    Rede,faixas = Configuration_Model(seed)

    Suscetiveis::UInt16 = S0;
    Expostos::UInt16 = E0;
    Assintomaticos::UInt16 = 0;
    Sintomaticos::UInt16 = 0;
    Hospitalizados::UInt16 = 0;
    Recuperados::UInt16 = 0;
    Mortos::UInt16 = 0;
    Nhospitalizados= 0
    N::UInt16 = length(faixas)

    prob_estagio = rand(N)
    estagio = zeros(UInt8,N)
    vacinado = zeros(UInt8,N)


    sintomatico = @SVector [29.1/100, 37.4/100, 41.68/100, 39.4/100, 31.3/100]

    hospitalizacao = @SVector [1.04/100, 1.33/100, 1.38/100, 7.6/100, 24/100]

    morte = @SVector [0.408/100, 1.04/100, 3.89/100, 9.98/100, 17.5/100]

    while E0 != 0
        @inbounds r = rand()
        N0 = Int(ceil(r*N))
        if estagio[N0] ==0
            estagio[N0] = 1
            E0 -= 1
        end
    end
    s = 1
    t::Float64 = 0.5
    ano = 0
    while true

        acumulada = cumsum([calc_estagio(site,Rede,estagio,vacinado,faixas,prob_estagio,morte,hospitalizacao) for site in 1:N])
        rate = acumulada[end]

        if rate == 0
            break
        end
        tempo = exponentialRand(rate)
        if tempo == 0
            break
        end
        ano += tempo
        Delta = rand()*rate
        rate = 0
        if ano >= 3*365
            break
        end

        site = searchsortedfirst(acumulada, Delta)

        if estagio[site] == 0 # Suscetiveis

            estagio[site] = 1
            Expostos += 1
            Suscetiveis -= 1

        elseif estagio[site] == 1 # Expostos
    
            sintoma = sintomatico[faixas[site]]*( vacinado[site] == 1 ? 0.34693877551 : 1)
            if prob_estagio[site] < sintoma
                estagio[site] = 3
                Sintomaticos += 1
            else
                estagio[site] = 2
                Assintomaticos += 1
            end
            Expostos -= 1
    
        elseif estagio[site] == 2 #Assintomáticos
    
            estagio[site] = 5
            Recuperados += 1
            Assintomaticos -= 1

    
        elseif estagio[site] == 3 # Sintomaticos
    
            if prob_estagio[site] < hospitalizacao[faixas[site]]
                estagio[site] = 4
                Hospitalizados += 1
                Nhospitalizados += 1
            else
                estagio[site] = 5
                Recuperados += 1
            end
            Sintomaticos -= 1
    
        elseif estagio[site] == 4
            if prob_estagio[site] < morte[faixas[site]]
                estagio[site] = 6
                Mortos += 1
            else
                estagio[site] = 5
                Recuperados += 1
            end
            Hospitalizados -= 1
    
        elseif estagio[site] == 5
            estagio[site] = 0
            Suscetiveis += 1
            Recuperados -= 1
        end

        prob_estagio[site] = rand()

        if ano > t 
            s += 1
            t += 0.5
        end
        grau = sum([length(neighbors(Rede,i)) for i in 1:N if estagio[i] == 1 || estagio[i] == 5 ])
        infect_time[s,1] += Suscetiveis;
        infect_time[s,2] += Expostos;
        infect_time[s,3] += Assintomaticos;
        infect_time[s,4] += Sintomaticos;
        infect_time[s,5] += Hospitalizados;
        infect_time[s,6] += Recuperados;
        infect_time[s,7] += Mortos;
        infect_time[s,8] += grau/(Suscetiveis+Recuperados)
        infect_time[s,9] += infect_time[s,8]*infect_time[s,8]
        quant[s] += 1

    end
    return [Nhospitalizados,Mortos,Nhospitalizados*Nhospitalizados,Mortos*Mortos],infect_time
end

function generate_infect(redes)

    N = length(readdlm("./input/faixas.txt",' '))
    tempo = 365*3*2;
    infect_time = zeros((tempo,9))
    quant = zeros(tempo)
    final = zeros(4)

    for i in 1:redes
        final_,infect_time = infect(i,N-1,1,infect_time,quant)
        final += final_
    end
    infect_time = [infect_time[i,:]/quant[i] for i in 1:tempo]
    #infect_time[:,9] = sqrt(infect_time[:,9] - infect_time[:,8]*infect_time[:,8])
    final = final./redes
    final[4] = final[4] - final[2]*final[2]
    final[3] = final[3] - final[1]*final[1]
end