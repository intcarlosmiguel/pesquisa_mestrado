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
    Q3
    Q5a_cQ4
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
    profissao = adultos['situation professionnelle '][int(i)]
    prof.append(profissao)
    if(profissao!=7):
        if(profissao>7):
            return prof
        else:
            a1 = float(adultos['Dans quel secteur d\'activité travaillez-vous'][int(i)])
            if(math.isnan(a1)):
                prof.append(0)
            else:
                prof.append(int(a1))
            a2 = float(adultos['profession qui entraîne beaucoup de contacts ?'][int(i)])
            if(math.isnan(a2)):
                prof.append(0)
            else:
                prof.append(int(a2))
            if(prof[-1]==1):
                prof.append([adultos['nombre de contacts en milieu pro'][int(i)],adultos["tranche(s) d'âge des contacts en milieu pro"][int(i)],adultos['A ou non plus de 20 contacts professionnels'][int(i)]])
            return prof
    else:
        a1 = adultos['Nombre d\'étudiants dans la classe'][int(i)]
        a2 = adultos['Etudiant qui mange  à la cantine ?'][int(i)]
        prof.append(a1)
        prof.append(a2)
        return prof

def create_adult(df,a,writer):
    adulto = {}
    adulto['Idade'] = df["Age du sujet de l'enquête"].values
    adulto['Sexo'] = df["Sexe du sujet de l'enquête"].values
    adulto['Casa'] = df["AGE PERSONNE"].values
    adulto["Escolaridade"] = df["diplôme le plus élevé"].values
    adulto["Locomocao_Semana"] = df["mode(s) deplacement semaine"].values
    adulto["Locomocao_Final_De_Semana"] = df["mode(s) deplacement we"].values
    adulto["Profissional"] = a
    adulto["Onda"] = df["vague"].values
    adulto["id"] = [i for i in range(df.index[-1]+1)]
    pd.DataFrame.from_dict(adulto).to_excel(writer, sheet_name='Adultos')

def create_children(df,a,writer):
    adulto = {}
    adulto['Idade'] = df["Age du sujet de l'enquête"].values
    adulto['Sexo'] = df["Sexe du sujet de l'enquête"].values
    adulto['Casa'] = df["AGE PERSONNE"].values
    adulto["Escolaridade"] = df["diplôme le plus élevé"].values
    adulto["Locomocao_Semana"] = df["mode(s) deplacement semaine"].values
    adulto["Locomocao_Final_De_Semana"] = df["mode(s) deplacement we"].values
    adulto["Escola"] = a
    adulto["Onda"] = df["vague"].values
    adulto["id"] = ['c'+str(i) for i in range(df.index[-1]+1)]
    pd.DataFrame.from_dict(adulto).to_excel(writer, sheet_name='Crianças')

def school(i,df):
    escola = df["L'enfant est-t-il scolarisé ?"][i]
    if(math.isnan(escola)):
        es = []
        es.append(2)
        return es
    else:
        es = [escola]
        es1 = []
        if(escola == 1):
            a1 = df["Nbr d'enfants dans sa classe"][i]
            a2 = df["enfant qui mange  à la cantine ?"][i]
            a3 = df["va au centre aéré ?"][i]
            es1.append(a1)
            es1.append(a2)
            es1.append(a3)
            if(a3==1):
                a4 = df["va au centre aéré durant ecole ?"][i]
                a5 = df["va au centre aéré durant vacances ?"][i]
                es1.append([a4,a5])
            es.append(es1)
            es.append([])
            return es
        if(escola == 2):
            a1 = df["enfant  gardé maison/en famille ?"][i]
            a2 = df["enfant gardé chez assist. Mat.\ nassistante maternelle ?"][i]
            es1.append(a1)
            es1.append(a2)
            if(a2== 1):
                a3 = df["nb enfants chez assistante"][i]
                a4 = df["l'assistante accueille des enft scolarisés ?"][i]
                es1.append([a3,a4])
            a1 = df["gardé en crèche ?"][i]
            es1.append(a1)
            if(a1==1):
                a2 = df["combien d'enfants dans la creche"][i]
                es1.append([a2])
            a3 = df["fréquence halte-garderie"][i]
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
    hist1 = locomotion(df,'Q5a_cQ4')
    hist2 = locomotion(df,'Q5b_cQ4b')

    c = ['indianred','steelblue','springgreen']
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

def stacked_bar(df,size1,size2,axs,name,number1,number2):
    font = {'family': 'monospace',
            'size': 14,
            }
    hist = np.zeros((size2,size1))
    for i in df:
        a = [(j[number1],j[number2]) for j in df[i]]
        for j in a:
            #print(a)
            if((not math.isnan(j[0])) and (not math.isnan(j[1]))):
                hist[int(j[0])-1][int(j[1])-1] += 1
    total = hist[0] + hist[1]
    hist = hist/total
    hist = hist.T
    cont = 0
    print(hist)
    for dados,nomes in zip(hist,name):
        axs.barh(nomes,dados[0] ,color = 'steelblue' ,label = 'Contato Físico' if cont == 0 else "")
        axs.barh(nomes,dados[1],left=dados[0],color = 'darkslategray',label = 'Sem Contato Físico'if cont == 0 else "" )
        cont += 1
    axs.legend()
    axs.set_ylabel('Duração',fontdict=font)
    return axs

def multiple_stacked_bar(df):

    fig, axs = plt.subplots(3,figsize=(12,10))
    fig.tight_layout(pad=5.0)
    axs[0] = stacked_bar(
        df,
        5,
        2,
        axs[0],
        ['< 5 minutos','5min -\n 15 min','15min - 1h','1h - 4h','> 4h'],
        6,
        -1
    )
    axs[1] = stacked_bar(
        df,
        5,
        2,
        axs[1],
        ['Diariamente','Semanalmente','Mensalmente','Anualmente','Primeira Vez'],
        6,
        5
    )
    axs[2] = stacked_bar(
        df,
        5,
        5,
        axs[2],
        ['Diariamente','Semanalmente','Mensalmente','Anualmente','Primeira Vez'],
        -1,
        5
    )
    #plt.legend()
    #axs[0].set_xlabel('Proporção de Conexões',fontdict=font)
    #plt.show()