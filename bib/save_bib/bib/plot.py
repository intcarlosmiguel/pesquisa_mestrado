import matplotlib.pyplot as plt
import numpy as np
import math
import plotly.graph_objs as go
from bib.cleaning import *
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib as mpl
from scipy.stats import chi2_contingency,ks_2samp

import plotly.express as px

import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots

font = {'family': 'monospace',
            'size': 16,
            }

def plot_mes(df,save):
    dic = {}

    for i,j in zip(df['dia/mês'].values,df['Número de Contatos'].values):
        if(not math.isnan(i)):
            r = str(int(i))
            dia = int(r[:-2])
            mês = int(r[-2:])
            if(mês in dic.keys()):
                dic[mês][dia-1] += j
            else:
                dic[mês] = [0]*31
                dic[mês][dia-1] += j
    fig, axs = plt.subplots(len(dic.keys()), figsize=(16, 16), )

    meses = [str(i) for i in list(dic.keys())]
    meses.sort()
    name =['Fevereiro','Março','Abril','Maio']

    font = {'family': 'monospace',
            'size': 15,
            }

    titulo = {'family': 'monospace',
            }

    for mês,count in zip( meses ,range((len(dic.keys())) )):
        a = range(1,32)
        x = dic[int(mês)]
        axs[count].bar(a,x,color = 'salmon')
        axs[count].set_ylabel("Número de Contatos",fontdict=font)
        axs[count].set_xlabel(f"{name[count]}",fontdict= font)
    fig.suptitle('Número de Contatos por Mês', fontdict=titulo,size = 18)
    plt.savefig(f'./img/{save}.jpg', dpi=300)
    plt.show()

def plot_contatos_idade(contatos,data,save,titulo):

    id, quantidade = np.unique(contatos['id'],return_counts=True)
    id = id.astype(int)
    df = {
        'id': id,
        'quantidade': quantidade
    }
    df = pd.DataFrame(df)
    df = pd.merge(data,df,how='outer',on='id')
    media_por_grupo = pd.DataFrame(df.groupby("Idade")["quantidade"].mean())
    
    plt.figure(figsize=(12,6))
    plt.bar(media_por_grupo.index,media_por_grupo['quantidade'],color = 'darkcyan')

    font = {'family': 'monospace',
            'size': 16,
            }

    plt.ylabel("Média de Contatos",fontdict=font)
    plt.xlabel("Idade",fontdict=font)
    plt.title(titulo,fontdict=font)
    plt.grid(True)

    plt.savefig(f'./img/{save}.jpg', dpi=300)
    plt.show()

def locomotion_hist(df,titulos,save,tam = 0.4):
    font = {'family': 'monospace',
            'size': 14,
            }
    hist1 = locomotion(df,'Q5a_cQ4a')
    hist2 = locomotion(df,'Q5b_cQ4b')

    c = ['indianred','steelblue','mediumseagreen']
    legenda = ['Carro Particular','Transporte Público','A pé']
    X_axis = np.arange(2)
    X = np.arange(-tam/2,tam,tam/2)[:-1]
    plt.figure(figsize=(6,10))

    for i in range(len(X)):
        plt.bar(X[i],hist1[i],tam,color = c[i],label = legenda[i])
        plt.bar(X[i]+3*tam,hist2[i],tam,color = c[i])
    
    plt.xticks([0,1], ['Semana','Finais de Semana\nou\nFeriados'],fontdict=font)
    plt.xlabel('Dias',fontdict=font)
    plt.ylabel("Número de Adultos",fontdict=font)
    plt.legend()
    plt.title(titulos,fontdict=font)
    plt.savefig(f'./img/{save}.jpg', dpi=300)
    plt.show()

def stacked_bar(df,size1,size2,axs,name,number1,number2,label,color,ylabel):

    font = {'family': 'monospace',
            'size': 14,
            }
    
    hist = np.zeros((size2,size1))

    for i in df:
        a = [(j[number1],j[number2]) for j in df[i]]
        for j in a:
            if((not math.isnan(j[0])) and (not math.isnan(j[1]))):
                hist[int(j[0])-1][int(j[1])-1] += 1
    
    total = np.sum(hist,axis = 0)
    hist = hist/total
    hist = hist.T
    cont = 0
    left = 0

    for dados,nomes in zip(hist,name):
        for i in range(len(dados)):
            axs.barh(nomes, dados[i], left = left ,color =  color[i],label = label[i] if cont == 0 else "")
            left += dados[i]
        cont += 1
        left = 0
    
    axs.legend()
    axs.set_ylabel(ylabel,fontdict=font)
    return axs

def multiple_stacked_bar(df,titulo,save):

    fig, axs = plt.subplots(3,figsize=(12,18))
    fig.tight_layout(pad=5.0)
    axs[0] = stacked_bar(
        df,
        5,
        2,
        axs[0],
        ['< 5 minutos','5min -\n 15 min','15min - 1h','1h - 4h','> 4h'],
        -2,
        -1,
        ['Contato Físico','Sem Conato Físico'],
        ['steelblue','darkslategray'],
        'Duração'
    )
    axs[1] = stacked_bar(
        df,
        5,
        2,
        axs[1],
        ['Diariamente','Semanalmente','Mensalmente','Anualmente','Primeira Vez'],
        -2,
        -3,
        ['Contato Físico','Sem Conato Físico'],
        ['steelblue','darkslategray'],
        'Frequência'
    )
    axs[2] = stacked_bar(
        df,
        5,
        5,
        axs[2],
        ['Diariamente','Semanalmente','Mensalmente','Anualmente','Primeira Vez'],
        -1,
        -3,
        ['< 5 minutos','5min -\n 15 min','15min - 1h','1h - 4h','> 4h'],
        ['darkslategray','teal','steelblue','dodgerblue','deepskyblue'][::-1],
        'Frequência'
    )
    fig.suptitle(titulo, fontsize=16)
    plt.savefig(f'./img/{save}.jpg', dpi=300)

def conncection_idade(contatos,data):

    total = pd.merge(data,contatos,how='outer',on='id')
    media = total.groupby("Idade")["idade"].mean().reset_index()
    mediana = total.groupby("Idade")["idade"].median().reset_index()

    plt.figure(figsize=(12,6),dpi = 500)
    plt.bar(media['Idade'],media['idade'],color = 'darkslategray',label = 'Média')
    plt.bar(mediana['Idade'],mediana['idade'],color = 'seagreen', label = 'Mediana')
    plt.grid(True)
    plt.legend()
    plt.xticks(np.arange(0, np.max(np.array(media.index)), 10))
    plt.xlabel('Idade dos Entrevistado',fontdict=font)
    plt.ylabel('Idade dos Contatos',fontdict=font)
    plt.savefig('./img/idade_idade.jpg')
    plt.show()

def heat_conncetion(data,contatos,contatos02):
    faixas = [
        [0,20],
        [20,30],
        [30,50],
        [50,70],
        [70,100000]
    ]
    data['Faixas'] = [check_faixa(i,faixas) for i in data['Idade'] ]
    contatos['faixas'] = [check_faixa(i,faixas) for i in contatos['idade'] ]
    contatos02['faixas'] = [check_faixa(i,faixas) for i in contatos02['idade'] ]
    novo_df = pd.merge(data, contatos, on='id', how='outer')[['Faixas','faixas']]
    novo_df2 = pd.merge(data, contatos02, on='id', how='outer')[['Faixas','faixas']]
    A = np.array(pd.crosstab(novo_df['Faixas'], novo_df['faixas']))
    B = np.array(pd.crosstab(novo_df2['Faixas'], novo_df2['faixas']))
    A_ = np.copy(A)
    B_ = np.copy(B)
    np.fill_diagonal(A,0)
    np.fill_diagonal(B,0)
    A = np.ceil((A+A.T)/2)
    B = np.ceil((B+B.T)/2)
    np.fill_diagonal(A,np.diag(A_/2))
    np.fill_diagonal(B,np.diag(B_/2))
    C = (A+B)/2
    C = np.ceil(C)
    return C

def adj(df,contacts,contacts02,Nmortos):

    faixa = [
            [0,20],
            [20,30],
            [30,50],
            [50,70],
            [70,100000]
        ]

    A = np.zeros((len(faixa),len(faixa)))
    a = 0
    for i,f in zip(faixa,range(len(faixa))):

        ids = df[(df['Idade'] >= i[0]) & (df['Idade'] < i[1])]['id'].values
        ids = [str(j) for j in ids]
        idades = list(contacts[contacts['id'].isin(ids)]['idade']) + list(contacts02[contacts02['id'].isin(ids)]['idade'])
        idades = [i for i in idades if(not math.isnan(i))]
        idades = [check_faixa(i,faixa) for i in idades]
        idades = np.array(idades)

        for j in range(len(faixa)):
            a += len(idades[idades==j])
            A[f][j] += len(idades[idades==j])
            A[f][j] = A[f][j]/len(ids)
    
    figure, ax = plt.subplots(figsize = (16,16),dpi = 500)

    ax.matshow(A, interpolation ='nearest',cmap = 'copper_r')
    #figure.colorbar(axes)

    for (i, j), z in np.ndenumerate(A):
        if(z>5):
            ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center',color = 'white')
        else:
            ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')

    xaxis = np.arange(len(list(Nmortos.keys())))
    lista = list(Nmortos.keys())
    lista[-1] = "Idade => 70"

    #lista[0] = "Idade < 1"
    lista = [i.replace(" <= ",r'$\leq$').replace("=>",r'$\geq$') for i in lista]

    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(lista)
    ax.set_yticklabels(lista)

    plt.title("Conexão relativa entre as faixas de idade",fontdict= font)
    plt.savefig("./img/map.jpg")
    plt.show()

def generate_vacinado(
        plot = 0,
        erro = 0,
        N = 7189,
        c = 0.0,
        save = True
    ):
    infect_vacinado = [i for i in os.listdir(f'./C/output/vacina/{N}/nponderado/') if(f"{c}0" in i)]
    prob = [float(i.split("_")[-1][:4]) for i in infect_vacinado]
    nomes = [i.split("_")[-2] for i in infect_vacinado]
    vac = []
    #cores = ["darkred","lightseagreen","darkolivegreen",'red','darkorange','royalblue','navy','purple','darkseagreen']
    #cores = ["#"+i for i in cores]
    for i in infect_vacinado:
        x = np.loadtxt(f'./C/output/vacina/{N}/nponderado/{i}')
        x = np.array([i for i in x if(sum(np.isnan(i)) == 0)]).T
        vac.append(x)
    integral = []
    fig = go.Figure()
    
    arquivo = {
        "kshell":"CK",
        "eigenvector":"CA",
        "clustering":"CC",
        "idade":"CI",
        "grau":"CG",
        "harmonic":"CH",
        "random":"CR",
        "close":"CP",
        "eccentricity":"CE",
        "betwenness":"CB",
        'graumorte': "CM",
        "probhosp":"PH",
        "probmorte":"PM",
        "probhospassin":"PHA",
        "probmortepassin":"PMA"
    }
    op = {
        "kshell":0.3,
        "eigenvector":0.3,
        "clustering":1,
        "idade":1,
        "grau":1,
        "harmonic":0.3,
        "random":1,
        "close":0.3,
        "eccentricity":1,
        "betweenness":0.3,
        #'graumorte': "CM"
    }

    df = {}
    M = []
    ylabel = []
    for infect,file,p in zip(vac,infect_vacinado,prob):
        integral.append(np.dot(infect[1 if(plot ==0) else 2], 0.01*np.ones(len(infect[1 if(plot ==0) else 2]))))
        y = infect[1 if(plot ==0) else 2]
        fig.add_trace(
            go.Scatter(
                x=infect[0][np.argsort(infect[0])], 
                y=y[np.argsort(infect[0])], 
                mode='markers', 
                
                name=arquivo[file.split("_")[0]],
                marker=dict(
                    size = 7,
                    #color = cores[file.split("_")[2]]
                ),
                line=dict(width=7),
                error_y=dict(
                    type='data',  # Ajuste para 'percent' se desejar barras de erro em percentagem
                    array=infect[3 if(plot ==0) else 4]/np.sqrt(200),
                    visible = erro == 1
                ),
                #opacity= op[file.split("_")[2]],
            )
        )
        #M.append(y[np.argsort(infect[0])][np.arange(1,101)%10 == 0])
        ylabel.append(arquivo[file.split("_")[0]])
        df[arquivo[file.split("_")[0]]] = infect[1 if(plot ==0) else 2]
    
    tit = 'Hospitalizados' if(plot ==0) else 'Mortos'
    fig.update_layout(
        width=800,  # Largura do gráfico em pixels
        height=800,  # Altura do gráfico em pixels
        xaxis=dict(title='Fração de Vacinados',tickfont=dict(size=25)),
        #xaxis=dict(),
        yaxis=dict(title='Fração de Hospitalizados'if(plot ==0) else 'Fração de Mortos',tickfont=dict(size=25)),
        template = "seaborn",
        #paper_bgcolor='rgba(0,0,0,0)',
        font=dict(
            #family="Courier New, monospace",
            size=30,
            #color="RebeccaPurple"
        ),
    )
    s = 20
    fig.update_layout(margin=dict(l=s, r=s, t=s, b=s))
    fig.show()
    
    
def generate_3D_graph(plot = 0,color = 'orange',cmap = 'bone_r'):

    infect_shape = [i for i in os.listdir("./C/output/infect/") if(('shape' in i) & ('txt' in i))]
    y = [
        "Eficácia contra Sintomáticos",
        "Eficácia para Espalhamento",
        "Tempo para se recuperar"
        ]
    for shape,y_ in zip(infect_shape,y):
        shape = np.loadtxt(f'./C/output/infect/{shape}').T
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(13.5, 5),dpi = 500, gridspec_kw={'width_ratios': [1, 0.7]})

        #fig = plt.figure(figsize=(8,8),dpi = 100)
        ax1.set_axis_off()
        ax1 = fig.add_subplot(121, projection='3d')
        # Criação do scatter plot em 3D
        ax1.scatter(shape[0], shape[1], shape[2 if(plot == 0) else 3], marker='o',s=10,c = color)

        # Definição dos rótulos dos eixos
        ax1.set_xlabel('Fração de Vacinados')
        ax1.set_ylabel(y_)
        ax1.set_zlabel('# de Hospitalizados' if(plot == 0) else '# de Mortos')
        #ax1.set_box_aspect([1, 1, 1])
        # Exibição do gráfico
        Nx= len(np.unique(shape[1]))
        Ny= len(np.unique(shape[0]))
        
        S = shape[2 if(plot == 0) else 3].reshape(101,-1)
        #ax2 = fig.add_subplot(122)
        #fig, ax = plt.subplots(figsize = (10,6))
        im = ax2.imshow(S, cmap=cmap)
        ax2.set_xticks(np.arange(Nx))
        ax2.set_xticklabels(['']*Nx)
        ax2.set_yticks(np.arange(Ny))
        ax2.set_yticklabels(['']*Ny)
        ax2.set_ylabel('Fração de Vacinados')
        ax2.set_xlabel('1 - Eficácia da Vacina')
        cbar = ax2.figure.colorbar(im, ax=ax2)

        pos = ax1.get_position()
        ax1.set_position([pos.x0, pos.y0 + 0.02, pos.width, pos.height])
        fig.suptitle('Gráfico 3D de Hospitalizados'if(plot == 0) else 'Gráfico 3D de Mortos')
        plt.savefig(f'./img/infect/shape_hospitalizados_{y_}.png' if(plot == 0) else f'./img/infect/shape_mortos_{y_}.png')
        plt.show()

def infect_distribution():
    vac = np.loadtxt("./C/output/infect/infect_0.50.txt").T

    plt.figure(figsize=(8,5),dpi=500)
    plt.errorbar(np.arange(3*365*2)/2,vac[-2],yerr = vac[-1], markersize=5, capsize=4,c = 'royalblue',ecolor = 'gray',marker = 'D')
    #
    #plt.legend()
    plt.grid()
    plt.ylabel('Grau Médio')
    plt.xlabel('Tempo')
    #plt.ylim(0.01,0.1)
    plt.xlim(0,60)
    plt.savefig('./img/infect/grau_tempo.jpg')
    plt.show()

def time():
    
    infect_time = [i for i in os.listdir("./C/time")]
    time = np.arange(3*365*2)/2

    mpl.rcParams['font.size'] = 16
    
    plt.figure(figsize=(12,8),dpi = 500)
    plt.style.use('ggplot')
    infect_without = np.loadtxt("./C/time/"+'infect_without.txt').T
    plt.scatter(time,infect_without[0],label = "Suscetível",c = 'seagreen')
    plt.scatter(time,infect_without[1],label = "Expostos",c = 'steelblue')
    plt.scatter(time,infect_without[5],label = "Recuperados",c = 'tomato')
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("Número de Sítios")
    #plt.ylim(0,750)
    plt.xlim(0,60)
    plt.savefig("./img/infect/primeiros_60_dias.png")
    plt.show()


    mpl.rcParams['font.size'] = 18
    fig, axs = plt.subplots(1, 2, figsize=(18, 8),dpi =500)
    #plt.figure(figsize=(12,8),dpi = 500)
    legenda = ['VI',
               'VG',
               'VR',
               ''
               ]
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        print(file)
        if('without' in file):
            axs[0].plot(time,infect[4],label = "Sem vacina",linestyle="--")
            axs[1].plot(time,infect[6],label = "Sem vacina",linestyle="--")
        else:
            axs[0].scatter(time,infect[4],label = "Vacinação "+leg,s = 50, marker= 'D')
            axs[1].scatter(time,infect[6],label = "Vacinação "+leg,s = 50, marker= 'D')
            
    plt.grid()
    #axs[0].legend()
    axs[1].legend()
    axs[0].set_xlabel("Dias")
    axs[0].set_ylabel("H(t)")
    axs[1].set_xlabel("Dias")
    axs[1].set_ylabel("D(t)")
    axs[0].set_ylim(0,10)
    #axs[1].set_ylim(0,200)
    axs[0].set_xlim(60,600)
    #axs[1].set_xlim(60,600)
    fig.tight_layout()
    plt.savefig("./img/infect/double_infeccao_tempo.png")
    plt.show()

    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(12,8),dpi = 500)
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        if('without' in file):
            plt.plot(time,infect[0],label = "Sem vacina",linestyle="--")
        else:
            plt.scatter(time,infect[0],label = "Vacinação "+leg,s = 5, marker= 'D',alpha = 0.7)
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("Número de Suscetíveis")
    plt.ylim(0,750)
    plt.xlim(60,600)
    #plt.savefig("./img/infect/infeccao_tempo.png")
    plt.show()


    plt.figure(figsize=(12,8),dpi = 500)
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        if('without' in file):
            plt.plot(time,infect[1],label = "Sem vacina",linestyle="--")
        else:
            plt.scatter(time,infect[1],label = "Vacinação "+leg,s = 20, marker= 'D',alpha = 0.7)
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("E(t)")
    plt.ylim(0,200)
    plt.xlim(60,600)
    plt.savefig("./img/infect/infeccao_tempo.png")
    plt.show()

    plt.figure(figsize=(12,8),dpi = 500)
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        if('without' in file):
            plt.plot(time,infect[3],label = "Sem vacina",linestyle="--")
        else:
            plt.scatter(time,infect[3],label = "Vacinação "+leg,s = 5, marker= 'D',alpha = 0.7)
        
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("Número de Assintomáticos")
    plt.ylim(0,100)
    plt.xlim(60,600)
    #plt.savefig("./img/infect/infeccao_tempo.png")
    plt.show()

    mpl.rcParams['font.size'] = 14
    
    plt.figure(figsize=(12,8),dpi = 500)
    plt.style.use('ggplot')
    infect_without = np.loadtxt("./C/time/"+'infect_without.txt').T
    infect_random = np.loadtxt("./C/time/infect_aleatorio.txt").T
    plt.scatter(time,infect_without[4],label = "Hospitalizados sem Vacina",c = 'tomato')
    plt.scatter(time,infect_random[4],label = "Hospitalizados com vacinação aleatória",c = 'steelblue')
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("# de Hospitalizados")
    #plt.ylim(0,200)
    plt.savefig("./img/infect/hospitalizados_aleatorio.png")
    plt.xlim(0,200)
    plt.show()

    plt.figure(figsize=(12,8),dpi = 500)
    plt.style.use('ggplot')
    plt.plot(time,infect_without[6],label = "Hospitalizados sem Vacina",c = 'tomato')
    plt.plot(time,infect_random[6],label = "Hospitalizados com vacinação aleatória",c = 'steelblue')
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("# de Mortos")
    plt.ylim(0,10)
    plt.savefig("./img/infect/mortos_aleatorio.png")
    plt.xlim(0,200)
    plt.show()

    plt.figure(figsize=(12,8),dpi = 500)
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        if('without' in file):
            plt.plot(time,infect[6],label = "Sem vacina",linestyle="--")
        else:
            plt.scatter(time,infect[6],label = "Vacinação "+leg,s = 8, marker= 'D',alpha = 0.7)
        reg = LinearRegression()
        reg.fit(time[time > 100].reshape(-1, 1), infect[6][time > 100])
        y_pred = reg.predict(time[time > 100].reshape(-1, 1))
        r2 = r2_score(infect[6][time > 100].reshape(-1, 1), y_pred)
        print(leg,[reg.coef_]+[reg.intercept_] ,r2)
    plt.grid()
    plt.legend()
    plt.xlabel("Dias")
    plt.ylabel("# de Mortos")
    plt.savefig("./img/infect/mortos_tempo.png")
    plt.show()

    

    # Treine o modelo com os dados
    

    plt.figure(figsize=(12,8),dpi = 500)
    for file,leg in zip(infect_time,legenda):
        infect = np.loadtxt("./C/time/"+file).T
        if('without' in file):
            plt.plot(time,infect[4],label = "Sem vacina",linestyle="--")
        else:
            plt.scatter(time,infect[4],label = "Vacinação "+leg,s = 5, marker= 'D',alpha = 0.7)
    plt.grid()
    plt.legend()
    plt.ylim(0,10)
    plt.xlim(60,600)
    plt.xlabel("Dias")
    plt.ylabel("# de Hospitalizados")
    plt.savefig("./img/infect/hospitalizado_tempo.png")
    plt.show()

def plot_idades(a):
    idades = a.drop_duplicates(subset='id')
    idades = idades['Idade'].values
    idades = idades[idades < 200]
    x = np.arange(np.max(idades)+1)
    hist = np.zeros(len(x))
    for i in idades:
        hist[int(i)] += 1
    hist = hist/np.sum(hist)
    plt.scatter(x,hist)
    plt.title("Histograma de Idades")
    plt.grid()
    plt.show()

def heat_map(
        Matrix,
        xlabel,
        ylabel,
        name,
        xname = 'Faixa etária',
        yname = 'Faixa etária',
        round_ = 3,
        percentage = False,
    ):
    fig = go.Figure(data=go.Heatmap(z=Matrix,colorscale = px.colors.sequential.Cividis_r))
    fontsize = 15
    for i in range(Matrix.shape[0]):
        for j in range(Matrix.shape[1]):
            fig.add_annotation(
                text=str(np.round(Matrix[i][j]*100,round_)) + "%" if(percentage) else str(np.round(Matrix[i][j],round_)),
                x=j,
                y=i,
                xref='x',
                yref='y',
                showarrow=False,
                font=dict(color="black" if(Matrix[i][j] < 0.72*np.max(Matrix)) else 'white')
            )
    fig.update_layout(
        height = 700,
        width = 700,
        xaxis=dict(
            title=xname,
            tickvals = np.arange(Matrix.shape[1]),
            ticktext = xlabel,
            tickfont=dict(size=fontsize)
        ),
        yaxis=dict(
            title=yname,
            tickvals = np.arange(Matrix.shape[0]),
            ticktext = ylabel,
            tickfont=dict(size=fontsize)
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        font=dict(
            #family="Courier New, monospace",
            size=fontsize,
            #color="RebeccaPurple"
        ),
    )
    s = 20
    fig.update_layout(margin=dict(l=s, r=s, t=s, b=s))
    fig.show()
    if(len(name) != 0):
        fig.write_image(f"./img/{name}.png")


def creates_plot(data,coluna1,coluna2,tr1,tr2,fig,row,col,showlegend =True,cor = px.colors.qualitative.T10):
    df = data[[coluna1,coluna2]]
    df = df.dropna()
    phy = df[coluna2].unique()
    dados = {
        coluna2 : [],
    }
    
    for i in tr2:
        dados[i] = []

    for i in phy:
        x,y = np.unique(df[df[coluna2] == i][coluna1],return_counts=True)
        dados[coluna2].append(tr1[i])
        y = y/np.sum(y)
        for k in range(len(tr2)):
            dados[tr2[k]].append(y[k])
    
    dados = pd.DataFrame(dados)
    chi2_stat, p_val, dof, expected = chi2_contingency(dados.values[:,1:])
    # Exibe os resultados
    #print("Estatística qui-quadrado:", chi2_stat)
    #print("Valor p:", p_val)
    #print("Graus de liberdade:", dof)
    #print("Frequências esperadas:")
    #print(expected)
    for i in range(1, len(dados.columns)):
        fig.add_trace(
            go.Bar(
                x=dados[coluna2][np.argsort(dados.iloc[:, i])],
                y=dados.iloc[:, i][np.argsort(dados.iloc[:, i])],
                name=dados.columns[i],
                showlegend=showlegend,
                marker=dict(color=cor[i-1]),
                legendgroup = row+col*2,
            ),
            
            row = row,
            col = col,
        )

    return fig


def create_multi_bar(polymod):
    fig = make_subplots(rows=2, cols=2)
    tr1 = {
        1:'< 5 minutos',
        2:'Entre 5 e 15 minutos',
        3:'Entre 15 minutos e 1 hora',
        4:'Entre 1 e 4 horas',
        5:'> 4 horas'
    }
    tr2 = ['Com contato físico','Sem contato físico']
    fig = creates_plot(polymod,'phys_contact','duration_multi',tr1,tr2,fig,1,1)
    tr1 = {
        1:'Diário',
        2:'Semanalmente',
        3:'Mensalmente',
        4:'1 vez por ano',
        5:'Primeira vez'
    }
    tr2 = ['Com contato físico','Sem contato físico']
    fig = creates_plot(polymod,'phys_contact','frequency_multi',tr1,tr2,fig,1,2,False)

    tr1 = {
        1:'Diário',
        2:'Semanalmente',
        3:'Mensalmente',
        4:'1 vez por ano',
        5:'Primeira vez'
    }

    tr2 = [
        '< 5 minutos',
        'Entre 5 e 15 minutos',
        'Entre 15 minutos e 1 hora',
        'Entre 1 e 4 horas',
        '> 4 horas'
    ]

    fig = creates_plot(polymod,'duration_multi','frequency_multi',tr1,tr2,fig,2,1,True,px.colors.qualitative.Prism)

    cor = px.colors.qualitative.T10
    df = polymod[['phys_contact','cnt_home','cnt_work','cnt_school','cnt_transport','cnt_leisure','cnt_otherplace']].dropna()
    dados = {
        'locais' : [
            'Casa',
            'Trabalho',
            "Escola",
            "Transporte",
            "Lazer",
            'Outros Locais'
        ],
        'Com contato físico': [],
        "Sem contato físico": [],
    }
    for coluna in df.columns[1:]:
        x,y = np.unique(df['phys_contact'].values[df[coluna].astype(bool).values],return_counts=True)
        y = y/np.sum(y)
        dados['Com contato físico'].append(y[0])
        dados['Sem contato físico'].append(y[1])
    dados = pd.DataFrame(dados)
    for i in range(1, len(dados.columns)):
            
        fig.add_trace(
            go.Bar(
                x=dados['locais'][np.argsort(dados.iloc[:, i])],
                y=dados.iloc[:, i][np.argsort(dados.iloc[:, i])],
                name=dados.columns[i],
                showlegend=False,
                marker=dict(color=cor[i-1]),
                legendgroup = 2+2*2,
            ),
            
            row = 2,
            col = 2,
        )

    fig.update_layout(
        width = 1000,
        height = 800,
        barmode='stack',
        template = 'seaborn',
        #legend_tracegroupgap=390
    )

    fig.show()
    fig.write_image("./img/graficos.png")


def compara(modelo_k,polymod):
    contagem = pd.crosstab(polymod['id'], polymod['Contato_idadeFaixas']).values
    faixas = polymod.drop_duplicates(subset='id')['IdadeFaixas'].values

    fig = make_subplots(rows=5, cols=1, subplot_titles= [f'Faixa {i+1}' for i in range(5)])
    color = px.colors.qualitative.Dark24
    for i in range(5):
        k = np.sum(contagem,axis = 1)[faixas == i]
        model_k = modelo_k.T[0][modelo_k.T[1] == i]
        x,y = np.unique(k,return_counts=True)
        x_,y_ = np.unique(model_k,return_counts=True)
        y = y/np.sum(y)
        y_ = y_/np.sum(y_)
        
        fig.add_trace(go.Bar(
                x=x,
                y=y,
                showlegend = True if(i == 0) else False,
                #mode='markers', 
                marker=dict(color=color[0]),
                name = 'Dados',
            ),
            row=i+1,
            col=1,
        )

        fig.add_trace(go.Bar(
                x=x,
                y=y_,
                showlegend = True if(i == 0) else False,
                #mode='markers', 
                marker=dict(color=color[1]),
                name = 'Modelo',
            ),
            row=i+1,
            col=1,
            
        )
        fig.update_xaxes(
            tickfont=dict(size=15),
            row=i+1,
            col=1,
            range=[0, 40]
        )
        statistic, p_value = ks_2samp(k, model_k)

        # Interpretando o resultado
        alpha = 0.05  # Nível de significância
        if p_value > alpha:
            print(f"As amostras provavelmente vêm da mesma distribuição {round(p_value,3)}.")
        else:
            print("As amostras provavelmente vêm de distribuições diferentes.")
    fig.update_layout(
        width=600,  # Largura do gráfico em pixels
        height=1000,  # Altura do gráfico em pixels
        #xaxis=dict(range=[0, 40],tickfont=dict(size=15)),
        #xaxis=dict(),
        #yaxis=dict(title='Fração',tickfont=dict(size=15)),
        template = "seaborn",
        paper_bgcolor='rgba(0,0,0,0)',
        font=dict(
            size=15,
        ),
    )
    s = 22
    fig.update_layout(margin=dict(l=s, r=s, t=s, b=s))
    fig.show()
    fig.write_image("./img/comparacao.png")