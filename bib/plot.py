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
    A = (A+A.T)/2
    B = (B+B.T)/2
    np.fill_diagonal(A,np.diag(A_/2))
    np.fill_diagonal(B,np.diag(B_/2))
    C = (A+B)/2
    C = C.astype(int)
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


def generate_vacinado(plot = 0,erro = 0):
    dir = [i for i in os.listdir('./C/output/infect/')]
    infect_vacinado = [i for i in dir if(('vacinado' in i) and ('0.50' not in i))]

    vac = []

    for i in infect_vacinado:
        x = np.loadtxt(f'./C/output/infect/{i}')
        x = np.array([i for i in x if(sum(np.isnan(i)) == 0)]).T
        vac.append(x)
    plt.figure(figsize=(8,5),dpi=500)
    
    for infect,file in zip(vac,infect_vacinado):
        if(erro == 0):
            plt.scatter(infect[0],infect[1 if(plot ==0) else 2],label = file,s = 5,marker ='D')
        elif(erro == 1):
             if(infect.shape[0] == 5):
                plt.errorbar(infect[0],infect[1 if(plot ==0) else 2],yerr = infect[3 if(plot ==0) else 4] , markersize=5,errorevery = 1,elinewidth= 1 ,linewidth = 0, capsize=5,marker = 'D',label = file)
    plt.plot()
    plt.legend()
    plt.grid()
    plt.ylabel('# de Hospitalizados'if(plot ==0) else '# de Mortos')
    plt.xlabel('Fração de Vacinados')
    plt.savefig('./img/infect/hospitalizado_tempo.jpg' if(plot ==0) else './img/infect/hospitalizado_tempo.jpg')
    plt.show()

def generate_3D_graph(plot = 0,color = 'orange',cmap = 'bone_r'):

    infect_shape = [i for i in dir if('shape' in i)]
    shape = np.loadtxt(f'./C/output/infect/{infect_shape[0]}').T
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(13.5, 5),dpi = 500, gridspec_kw={'width_ratios': [1, 0.7]})

    #fig = plt.figure(figsize=(8,8),dpi = 100)
    ax1.set_axis_off()
    ax1 = fig.add_subplot(121, projection='3d')
    # Criação do scatter plot em 3D
    ax1.scatter(shape[0], shape[1], shape[2 if(plot == 0) else 3], marker='o',s=10,c = color)

    # Definição dos rótulos dos eixos
    ax1.set_xlabel('Fração de Vacinados')
    ax1.set_ylabel('1 - Eficácia da Vacina')
    ax1.set_zlabel('# de Hospitalizados' if(plot == 0) else '# de Mortos')
    #ax1.set_box_aspect([1, 1, 1])
    # Exibição do gráfico

    S = shape[2 if(plot == 0) else 3].reshape(101,-1)
    #ax2 = fig.add_subplot(122)
    #fig, ax = plt.subplots(figsize = (10,6))
    im = ax2.imshow(S, cmap=cmap)
    ax2.set_xticks(np.arange(101))
    ax2.set_xticklabels(['']*101)
    ax2.set_yticks(np.arange(101))
    ax2.set_yticklabels(['']*101)
    ax2.set_ylabel('Fração de Vacinados')
    ax2.set_xlabel('1 - Eficácia da Vacina')
    cbar = ax2.figure.colorbar(im, ax=ax2)

    pos = ax1.get_position()
    ax1.set_position([pos.x0, pos.y0 + 0.02, pos.width, pos.height])
    fig.suptitle('Gráfico 3D de Hospitalizados'if(plot == 0) else 'Gráfico 3D de Mortos')
    plt.savefig('./img/infect/shape_hospitalizados.png' if(plot == 0) else './img/infect/shape_mortos.png')
    plt.show()