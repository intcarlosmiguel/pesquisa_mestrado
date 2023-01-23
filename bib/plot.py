import matplotlib.pyplot as plt
import numpy as np
import math
from bib.cleaning import *

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

def plot_contatos_idade(contatos,df,save,titulo):

    idades = df["Q1_cQ1"].values
    contatos = contatos["Número de Contatos"].values
    hist = list(np.arange(np.max(idades)+1))
    value = np.zeros(len(hist))

    for idade,c in zip(idades,contatos):
        value[idade] += c
    
    plt.figure(figsize=(12,6))
    plt.bar(hist,value,color = 'darkcyan')

    font = {'family': 'monospace',
            'size': 16,
            }

    plt.ylabel("Número de Contatos",fontdict=font)
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

def conncection_idade(contatos,df):

    id_ = []
    a = 0
    c = 0
    for i in df['Q1_cQ1'].values:
        if(i>15):
            id_.append(str(a))
            a += 1
        else:
            id_.append('c' + str(c))
            c += 1
    
    id_ = np.array(id_)
    idades = df["Q1_cQ1"].values
    hist = list(np.arange(np.max(idades)+1))
    value = []
    for idade in hist:
        teste = []
        for i in id_[idades==idade]:
            if(i in contatos):
                teste.append([j[0] for j in contatos[i] if(not math.isnan(j[0]))])
        value.append([j for i in teste for j in i])
    
    media = [np.mean(i) for i in value]
    mediana = [np.median(i) for i in value]

    plt.figure(figsize=(12,6))
    plt.bar(hist,media,color = 'darkslategray',label = 'Média')
    plt.bar(hist,mediana,color = 'seagreen', label = 'Mediana')
    plt.grid(True)
    plt.legend()
    plt.xticks(np.arange(0, np.max(np.array(hist)), 10))
    plt.xlabel('Idade dos Entrevistado',fontdict=font)
    plt.ylabel('Idade dos Contatos',fontdict=font)
    plt.show()

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