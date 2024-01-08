import networkx as nx
import numpy as np
from IPython.display import clear_output
from itertools import combinations

def calculate(G):

    degree2 = np.array([G.degree[i] for i in G.nodes])
    media = np.mean(degree2)
    mediana = np.median(degree2)
    std = np.std(degree2)
    clustering = nx.average_clustering(G)
    path = nx.average_shortest_path_length(G)

    cl = np.array(list(nx.clustering(G).values()))
    r = np.corrcoef(degree2, cl)
    diametro = nx.diameter(G)
    
    return np.array([media,mediana,std,clustering,r[0][1],path,diametro])

def sorting(array,sorted_by):
    sort_ = np.argsort(sorted_by)[::-1]
    return array[sort_],sorted_by[sort_]

def conf_model_p(G,deg,p,site,n_existir):

    vizinhos = np.array(list(nx.neighbors(G,site)))
    vizinhos = vizinhos[deg[vizinhos] > 0]
    if(len(vizinhos)> 1):
        
        for i in vizinhos:

            shuff = vizinhos.copy()
            shuff = shuff[shuff != i]
            np.random.shuffle(shuff)

            for j in shuff:
                if(deg[i] == 0):
                    break

                rand = np.random.random(1)[0]
                if(rand <= p):
                    if((deg[j] == 0) or G.has_edge(i, j)):
                        continue
                    if(((i,j) in n_existir) or ((j,i) in n_existir)):
                        continue
                    deg[i] -= 1
                    deg[j] -= 1
                    G.add_edge(i,j)
                else:
                    if(((i,j) not in n_existir) or ((j,i) not in n_existir)):
                        n_existir.append((i,j))
    return deg,n_existir

def generate_config_model(T,degree,seed,p = 0):
    resultado = np.array([configuration_model(degree,i+seed,i,p) for i in range(T)])
    if(T!=1):
        clear_output(wait=True)
    return np.sum(resultado,axis = 0)/T
    

def visualize(degree,T):
    p = [1.0]
    zeed = [3500,1000,2000]
    metrics = [generate_config_model(T,degree,seed,probabilidade) for probabilidade,seed in zip(p,zeed)]
    for i in range(len(metrics)):
        metrics[i] = ['{0:.3g}'.format(i) for i in metrics[i]]
    for metrica,probabilidade in zip(metrics,p):
        print(probabilidade,*metrica, sep = "\t")


def configuration_model(degree,seed,number,p = 0):
    
    deg = list(degree.copy())
    np.random.seed(seed)

    deg.sort()
    deg = np.array(deg[::-1])
    n_existir = []
    clear_output(wait=True)
    G = nx.Graph()
    for i in range(len(deg)):
        shuff = np.arange(len(deg))
        shuff = shuff[deg != 0]
        shuff = shuff[shuff != i]
        np.random.shuffle(shuff)
        #print(f"{i}/{len(deg)}")
        for j in shuff:
            if(deg[i] == 0):
                deg,n_existir = conf_model_p(G,deg,p,i,n_existir)
                break
            if((deg[j] == 0) or G.has_edge(i, j)):
                continue
            deg[i] -= 1
            deg[j] -= 1
            G.add_edge(i,j)
    
    if(not nx.is_connected(G)):
        largest_cc = max(nx.connected_components(G), key=len)
        return configuration_model(degree,seed+1,number,p)
    print(number,p)
    return calculate(G)

def ER(degree,T):
    resultado = []

    def connected(degree,i):
        seed = 0
        G = nx.gnm_random_graph(len(degree), np.sum(degree)/2, seed=i)
        while(not nx.is_connected(G)):
            G = nx.gnm_random_graph(len(degree), np.sum(degree)/2, seed=i+seed)
            seed += 500
        if((i+1)%10 == 0):
            print(i+1)
        return calculate(G)

    resultado = np.array([connected(degree,i) for i in range(T)])

    resultado = np.array(resultado)
    resultado = np.sum(resultado,axis = 0)/T
    metrics = ['{0:.3g}'.format(i) for i in resultado]
    print(*metrics, sep = "\t")
