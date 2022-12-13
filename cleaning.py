import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def change_age(data):
    idade = [i for i in data.columns if('AGE' in i)]
    idades = [list(i) for i in data[idade].fillna(-1).values]
    idades = [i[:i.index(-1)] if(-1 in i) else i for i in idades]
    for i in idade:
        data.drop(i, inplace=True, axis=1)
    data['AGE PERSONNE'] = idades
    return data

def change_mode(data):
    modes = [i for i in data.columns if('mode(s)' in i)]
    for i in range(len(modes)):
        mode = list(data.columns).index(modes[i])
        mode = data.columns[mode:mode+3]
        modes_ = [list(i) for i in data[mode].fillna(-1).values]
        modes_ = [i[:i.index(-1)] if(-1 in i) else i for i in modes_]
        for j in mode:
            data.drop(j, inplace=True, axis=1)
        data[modes[i]] = modes_
    return data 

def change_traches(data):
    tranches = [i for i in data.columns if('tranche(s)' in i)][0]
    tranche = list(data.columns).index(tranches)
    tranche_ = data.columns[tranche:tranche+5]
    tranche = [list(i) for i in data[tranche_].fillna(0).values]
    for i in tranche_:
        data.drop(i, inplace=True, axis=1)
    data[tranches] = tranche
    return data

def change_column(data):
    a = '''onda
    pergunta de número
    CONTROL Num Questionio versus recrutar Tel
    Tipo de questionário
    Q1_cQ1
    Q2
    cQ2
    cQ1
    cQ6
    Q_
    Q_1
    Q4_cQ8
    Q6_cQ9
    Q7_cQ10
    cQ11
    cQ12
    cQ13
    cQ13a
    cQ13b
    cQ14
    cQ14a
    cQ15
    cQ16
    cQ17
    cQ18
    cQ18a
    cQ18b
    Q8
    Q8a
    Q8c
    Q9
    Q10
    Q3_cQ7
    Q5a_cQ4a
    Q5b_cQ4b
    Q8b
    '''
    colunas = a.split('\n')
    colunas = [i.replace('    ','') for i in colunas]
    data.columns = colunas[:-1]
    return data

def remove_adultos(df,aux):
    if(aux==0):
        remocao = [i for i in df.columns if i[0]=='c']
    else:
        remocao = [i for i in df.columns if ((i[0]=='Q') & ('_c' not in i))]
    for i in remocao:
        df = df.drop(i, axis=1)
    df = df.drop("CONTROL Num Questionio versus recrutar Tel",axis = 1)
    df = df.drop("Tipo de questionário",axis = 1)
    return df

def profissional(i,adultos):
    prof = []
    profissao = adultos['Q6_cQ9'][int(i)]
    prof.append(profissao)
    if(profissao!=7):
        if(profissao>7):
            return prof
        else:
            a1 = float(adultos['Q7_cQ10'][int(i)])
            if(math.isnan(a1)):
                prof.append(0)
            else:
                prof.append(int(a1))
            a2 = float(adultos['Q8'][int(i)])
            if(math.isnan(a2)):
                prof.append(0)
            else:
                prof.append(int(a2))
            if(prof[-1]==1):
                prof.append([adultos['Q8a'][int(i)],adultos["Q8b"][int(i)],adultos['Q8c'][int(i)]])
            return prof
    else:
        a1 = adultos['Q9'][int(i)]
        a2 = adultos['Q10'][int(i)]
        prof.append(a1)
        prof.append(a2)
        return prof

def create_adult(df,a,writer):
    adulto = {}
    adulto['Idade'] = df["Q1_cQ1"].values
    adulto['Sexo'] = df["Q2"].values
    adulto['Casa'] = df["Q3_cQ7"].values
    adulto["Escolaridade"] = df["Q4_cQ8"].values
    adulto["Locomocao_Semana"] = df["Q5a_cQ4a"].values
    adulto["Locomocao_Final_De_Semana"] = df["Q5b_cQ4b"].values
    adulto["Profissional"] = a
    adulto["Onda"] = df["onda"].values
    adulto["id"] = [i for i in range(df.index[-1]+1)]
    pd.DataFrame.from_dict(adulto).to_excel(writer, sheet_name='Adultos',index = False)

def create_children(df,a,writer):
    adulto = {}
    adulto['Idade'] = df["Q1_cQ1"].values
    adulto['Sexo'] = df["cQ2"].values
    adulto['Casa'] = df["Q3_cQ7"].values
    adulto["Escolaridade"] = df["Q4_cQ8"].values
    adulto["Locomocao_Semana"] = df["Q5a_cQ4a"].values
    adulto["Locomocao_Final_De_Semana"] = df["Q5b_cQ4b"].values
    adulto["Escola"] = a
    adulto["Onda"] = df["onda"].values
    adulto["id"] = ['c'+str(i) for i in range(df.index[-1]+1)]
    pd.DataFrame.from_dict(adulto).to_excel(writer, sheet_name='Crianças',index= False)

def school(i,df):
    escola = df["cQ11"][i]
    if(math.isnan(escola)):
        es = []
        es.append(2)
        return es
    else:
        es = [escola]
        es1 = []
        if(escola == 1): # Se estiver na escola
            a1 = df["cQ16"][i]
            a2 = df["cQ17"][i]
            a3 = df["cQ18"][i]
            es1.append(a1)
            es1.append(a2)
            es1.append(a3)
            if(a3==1):
                a4 = df["cQ18a"][i]
                a5 = df["cQ18b"][i]
                es1.append([a4,a5])
            es.append(es1)
            es.append([])
            return es
        if(escola == 2):
            a1 = df["cQ12"][i]
            a2 = df["cQ13"][i]
            es1.append(a1)
            es1.append(a2)
            if(a2== 1):
                a3 = df["cQ13a"][i]
                a4 = df["cQ13b"][i]
                es1.append([a3,a4])
            a1 = df["cQ14"][i]
            es1.append(a1)
            if(a1==1):
                a2 = df["cQ14a"][i]
                es1.append([a2])
            a3 = df["cQ15"][i]
            es1.append(a3)
            es.append([])
            es.append(es1)
        return es

def get_number(vetor,value):
    number = np.arange(0,len(vetor),1)
    return number[vetor==value][0]

def get_contacts(k,i,contatos):
    contacts = []
    for j in range(i,i+4):
        contacts.append(contatos.iloc[k,j])
    contacts.append(list(contatos.iloc[0,i+4:i+11].values))
    for j in range(i+11,i+14):
        contacts.append(contatos.iloc[k,j])
    return contacts

def get_json(contatos,df):

    number = [get_number(contatos.columns,i) for i in contatos.columns if 'age min' in i]
    age = ['adulto' if i==1 else 'crianca' for i in df['Tipo de questionário']]
    cont1 = 0
    cont2 = 0
    age2 = []
    for a in age:
        if(a=='crianca'):
            age2.append('c'+str(cont1))
            cont1 += 1
        else:
            age2.append(str(cont2))
            cont2 += 1

    contacts = {}

    for z,j in zip(age2,range(contatos.shape[0])):
        contacts[str(z)] = [get_contacts(j,i,contatos) for i in number if(not math.isnan(contatos.iloc[j,i]))]
    return contacts

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

def plot_contatos_idade(df,criancas,adultos,save,titulo):
    idade_contatos = []
    for i in df:
        if(i[0] == 'c'):
            idade_contatos.append([criancas[criancas['index'] == i]['Q1_cQ1'].values[0],len(df[i])])
        else:
            idade_contatos.append([adultos[adultos['index'] == int(i)]['Q1_cQ1'].values[0],len(df[i])])

    hist = []
    value = []
    for i in idade_contatos:
        if(i[0] in hist):
            value[hist.index(i[0])] += i[1]
        else:
            hist.append(i[0])
            value.append(i[1])
    
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

def locomotion(df,coluna):
    hist = [0]*3
    for i in df[coluna].values:
        for j in range(len(i)):
            hist[j] += i[j]
    hist.sort()
    return np.array(hist)

def locomotion_adultos(df,titulos,save,tam = 0.4):
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
            #axs.barh(nomes,dados[1],left=dados[0],color = ,label = 'Sem Contato Físico'if cont == 0 else "" )
        cont += 1
        left = 0
    axs.legend()
    axs.set_ylabel(ylabel,fontdict=font)
    return axs

def multiple_stacked_bar(df,titulo):

    fig, axs = plt.subplots(3,figsize=(12,18))
    fig.tight_layout(pad=5.0)
    axs[0] = stacked_bar(
        df,
        5,
        2,
        axs[0],
        ['< 5 minutos','5min -\n 15 min','15min - 1h','1h - 4h','> 4h'],
        6,
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
        6,
        5,
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
        5,
        ['< 5 minutos','5min -\n 15 min','15min - 1h','1h - 4h','> 4h'],
        ['darkslategray','teal','steelblue','dodgerblue','deepskyblue'][::-1],
        'Frequência'
    )
    fig.suptitle(titulo, fontsize=16)
    #plt.legend()
    #axs[0].set_xlabel('Proporção de Conexões',fontdict=font)
    #plt.show()

def conncection_idade(df,criancas,adultos):

    idade = []
    idade_contatos = []
    hist = []
    value = []

    for i in df:
        idade_contatos.append([j[2] for j in df[i]])
        if(i[0] == 'c'):
            idade.append(criancas[criancas['index'] == i]['Q1_cQ1'].values[0])
        else:
            idade.append(adultos[adultos['index'] == int(i)]['Q1_cQ1'].values[0])
    

    for i,j in zip(idade,idade_contatos):
        if(i in hist):
            value[hist.index(i)].append(j)
        else:
            hist.append(i)
            value.append([j])

    for i in range(len(value)):
        value[i] = np.array([item for sublist in value[i] for item in sublist])
    
    media = [np.mean(i) for i in value]
    mediana = [np.median(i) for i in value]

    font = {'family': 'monospace',
            'size': 16,
            }
            
    plt.figure(figsize=(12,6))
    plt.bar(hist,media,color = 'darkslategray',label = 'Média')
    plt.bar(hist,mediana,color = 'seagreen', label = 'Mediana')
    plt.grid(True)
    plt.legend()
    plt.xticks(np.arange(0, np.max(np.array(hist)), 10))
    plt.xlabel('Idade dos Entrevistado',fontdict=font)
    plt.ylabel('Idade dos Contatos',fontdict=font)
    plt.show()


def check_faixa(idade,faixas):
    for faixa in range(len(faixas)):
        if(faixa == 0):
            if(idade < 1):
                return int(faixa)
        else:
            i = faixas[faixa]
            if((idade<=i[1]) and(idade>=i[0])):
                return int(faixa)

def calc_adj(A,adultos,faixa,contacts):
    for adulto in adultos['index']:
        idade = check_faixa(adultos['Q1_cQ1'][adultos['index'] == adulto].values,faixa)
        adulto = str(adulto)
        if(contacts[adulto] == []):
            continue
        idades = np.array([check_faixa(i[2],faixa) for i in contacts[adulto]])
        for i in idades:
            A[idade][i] += 1
    return A

def adj(adultos,criancas,contacts,contacts02,Nmortos):

    faixa = [
        [0,0],
        [1,5],
        [6,19],
    ]
    [faixa.append([i,i+9]) for i in range(20,100,10)]
    faixa[-1] = [90,100000]

    A = np.zeros((len(faixa),len(faixa)))
    A = calc_adj(A,adultos,faixa,contacts)
    A = calc_adj(A,adultos,faixa,contacts02)
    A = calc_adj(A,criancas,faixa,contacts)
    A = calc_adj(A,criancas,faixa,contacts02)

    for i in range(len(A)):
        if(np.sum(A[i])!= 0):
            A[i] = A[i]/np.sum(A[i])

    figure, ax = plt.subplots(figsize = (28,28),dpi = 300)

    axes = ax.matshow(A, interpolation ='nearest',cmap = 'copper_r')
    figure.colorbar(axes)

    for (i, j), z in np.ndenumerate(A):
        if(z>0.35):
            ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center',color = 'white')
        else:
            ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')

    xaxis = np.arange(len(list(Nmortos.keys())))
    lista = list(Nmortos.keys())
    lista[-1] = "Idade => 90"
    lista[0] = "1 > Idade"

    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(lista)
    ax.set_yticklabels(lista)

    plt.title("Conexão entre as faixas de idade")
    plt.savefig("./img/map.jpg")
    plt.show()

    return A