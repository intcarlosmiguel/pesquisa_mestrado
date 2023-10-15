import math
import numpy as np
import pandas as pd
import json
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def generate_contatos(contatos):
    colunas = [list(contatos.columns).index(i) for i in contatos.columns if('moyen' in i )]
    colunas.append(colunas[-1]+12)
    s = [list(contatos.columns).index(i) for i in contatos.columns if('Nombre de contacts' in i )][0]
    s = contatos.columns[s]
    numeros = contatos[s].values
    lista = []
    ids = []
    contato = {}
    c = []
    for i in range(1,len(colunas)):
        c += list(range(colunas[i-1],colunas[i]))
    df = contatos.iloc[:,c]
    df.insert(0, 'id', list(contatos.index))
    df = df.values
    dic = {}
    for row in df:
        dic[int(row[0])] = row[1:]
    # Reshape para a forma desejada
    #df = df.reshape((81320, 12))
    a = []
    for row in dic:
        dic[row] = dic[row].reshape(-1,12)
        d = pd.DataFrame(dic[row])
        d = d.dropna(subset=[0])
        d = d.fillna(-1).values
        a += np.hstack((np.ones(d.shape[0]).reshape(-1,1)*int(row),np.array(d))).tolist()
    a = np.array(a)
    
    """ for n,id in zip(numeros,contatos.index):
        if(n!=0):
            data = contatos.iloc[id,list(range(colunas[0],colunas[n]))]
            lista.append(np.reshape(data.values,(n,-1)))
            ids.append(id)
    
    lista = np.array([np.concatenate(([id],j)).astype(int) for i,id in zip(lista,ids) for j in i]) """

    contato['id'] = a.T[0].astype(int)
    contato['idade'] = a.T[1]
    contato['sexo'] = a.T[2]
    contato['pele'] = a.T[-2]
    contato['frequência'] = a.T[-3]
    contato['duração'] = a.T[-1]
    contato = pd.DataFrame.from_dict(contato)
    contato['local'] = list(a[:,3:10].astype(int))
    


    return contato

def generate_df():
    if('participantes.csv' not in os.listdir("./output/")):
        df = pd.read_excel('./input/RawData_ComesF.xlsx',index_col=False,skiprows=[0, 1])
        participantes = df.iloc[:,:54]
        contatos = df.iloc[:,54:]
        for i in contatos.columns:
            if(('min' in i) or ('max' in i)):
                del contatos[i]
        degree = contatos['Nombre de contacts saisis JOUR 1 + JOUR 2'].values
        np.savetxt('./output/degree.txt',np.ceil(degree[degree>0]/2),fmt='%d')
        x = list(contatos.columns).index('jour2mois2')
        contatos01 = contatos[contatos.columns[:x]]
        contatos02 = contatos[contatos.columns[x:]]
        contatos01 = generate_contatos(contatos01)
        contatos02 = generate_contatos(contatos02)
        data = participantes[[
            "Age du sujet de l'enquête",
            "Sexe du sujet de l'enquête",
            'Nombre de personnes au foyer',
            'situation professionnelle ',
            "Dans quel secteur d'activité travaillez-vous"
            ]
        ]
        data['id'] = list(df.index)
        colunas = ['Idade','Sexo','#Familiares','Profissão','Setor']
        for i,j in zip(data.columns,colunas):
            data = data.rename(columns={i:j})
        data['#Contatos01'] = contatos['Nombre de contacts saisis JOUR 1']
        data['#Contatos02'] = contatos['Nombre de contacts saisis JOUR 2']
        data['Dia01'] = contatos['jour1mois1']
        data['Dia02'] = contatos['jour2mois2']

        faixas = [
            [0,20],
            [20,30],
            [30,50],
            [50,70],
            [70,100000]
        ]
        data['Faixas'] = [check_faixa(i,faixas) for i in data['Idade'] ]
        contatos01['faixas'] = [check_faixa(i,faixas) for i in contatos01['idade'] ]
        contatos02['faixas'] = [check_faixa(i,faixas) for i in contatos02['idade'] ]

        pd.DataFrame.to_csv(data,'./output/participantes.csv',index=False)
        pd.DataFrame.to_csv(contatos01,f'./output/contatos_01.csv',index=False)
        pd.DataFrame.to_csv(contatos02,f'./output/contatos_02.csv',index=False)
        return data,contatos01,contatos02
    else:
        data = pd.read_csv('./output/participantes.csv')
        contatos01 = pd.read_csv('./output/contatos_01.csv')
        contatos02 = pd.read_csv('./output/contatos_02.csv')
        return data,contatos01,contatos02
    

def change_age(data):
    idade = [i for i in data.columns if('AGE' in i)]
    idades = [list(i) for i in data[idade].fillna(-1).values]
    idades = [i[:i.index(-1)] if(-1 in i) else i for i in idades]
    for i in idade:
        data.drop(i, inplace=True, axis=1)
    data['AGE PERSONNE'] = idades
    return data

def create_design(idades,contatos,contatos02):
    design1 = []
    design2 = []
    for i,j in zip(contatos02['dia/mês'].values,idades):
        if(not math.isnan(i)):
            mes = int(str(int(i))[-2:])
            if((mes==2)|(mes==3)):
                design1.append(j)
            if((mes==4) | (mes==5)):
                design2.append(j)
    for i,j in zip(contatos['dia/mês'].values,idades):
        if(not math.isnan(i)):
            mes = int(str(int(i))[-2:])
            if((mes==2)|(mes==3)):
                if(j not in design1):
                    design1.append(j)
            if((mes==4) | (mes==5)):
                if(j not in design2):
                    design2.append(j)
            #print(mes)
    np.savetxt('./dados/design1.txt',design1,fmt='%s')
    np.savetxt('./dados/design2.txt',design2,fmt='%s')

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

def locomotion(df,coluna):
    hist = [0]*3
    for i in df[coluna].values:
        for j in range(len(i)):
            hist[j] += i[j]
    hist.sort()
    return np.array(hist)

def check_faixa(idade,faixas):
    for faixa in range(len(faixas)):
        i = faixas[faixa]
        if((idade<i[1]) and(idade>=i[0])):
            return int(faixa)

def generate_mortalidade(df,name = 'mortalidade',vacina = 'nenhuma'):

    hist = [
        [0,20],
        [20,30],
        [30,50],
        [50,70],
        [70,100000]
    ]
    label = [f'{i[0]} <= Idade < {i[1]} ' for i in hist]

    # Número de mortos dividido por Número de Infectados
    Nmortos = {}
    for i,j in zip(label,hist):
        infectados = df[(df['NU_IDADE_N']>=j[0]) & (df['NU_IDADE_N']<=j[1])]
        mortos = 0
        if(vacina == 'nenhuma'):
            mortos = infectados[infectados['EVOLUCAO']==2]
        else:
            m = infectados[infectados['FAB_COV_1'].str.contains(vacina.upper())]
            mortos = m[m['EVOLUCAO']==2]
        Nmortos[i] = mortos.shape[0]/infectados.shape[0]
    
    with open(f'./dados/{name}.json', 'w') as f:
        json.dump(Nmortos, f)
    return Nmortos

def where(vetor,condition):
    a = np.arange(len(vetor))
    return a[condition]

def contacts_to_df(contacts):
    cont = [[i]*len(contacts[i]) for i in contacts]
    cont = [int(j) for i in cont for j in i]

    idades = [j[0]for i in contacts for j in contacts[i]]
    freq = [j[3]for i in contacts for j in contacts[i]]
    pele = [j[4]for i in contacts for j in contacts[i]]
    duracao = [j[5]for i in contacts for j in contacts[i]]

    contatos = {}
    contatos['id'] = cont

    contatos['idade'] = idades
    contatos['Frequência'] = freq
    contatos['Pele'] = pele
    contatos['Duração'] = duracao
    contatos = pd.DataFrame.from_dict(contatos)
    return contatos

def create_degree(data):
    degree01 = data['#Contatos01'].values
    degree02 = data['#Contatos02'].values
    np.savetxt('./dados/degree01.txt',degree01,fmt = "%d")
    np.savetxt('./dados/degree02.txt',degree02,fmt = "%d")

def LM(x,y,axs,color = 'reds'):
    X = x.reshape(-1,1)
    regressor = LinearRegression()
    regressor.fit(X, y)
    coeficientes = regressor.coef_[0]
    intercepto = regressor.intercept_

    # Calcular o coeficiente de determinação (R²)
    r2 = regressor.score(X, y)
    if(axs != -1):
        axs.plot(x,x*coeficientes + intercepto,c = color,label = 'R² :{:.3f}, {:.3f}x + {:.3f}'.format(r2,coeficientes,intercepto))
    return coeficientes,intercepto,r2

def degree_distribution(graus,faixas,tipo = 0,grafico = 0):
    c = ['red','blue','green','k','navy']
    M = np.zeros((5,5))
    A = np.zeros((5,5))
    for i in range(5):
        degree = graus[faixas == i,:].T
        rang = np.arange(np.max(degree)+1)
        if(grafico == 1):
            fig, axs = plt.subplots(5,figsize=(8, 20))
        k = 0
        for deg,color in zip(degree,c):
            hist = np.zeros(len(rang))
            for j in deg:
                hist[int(j)] += 1
            hist = hist/np.sum(hist)
            x = rang[hist > 0]
            hist = hist[hist > 0]
            
            hist = hist[x<10]
            x = x[x<10]

            if(tipo == 1):
                hist = hist[x > 2]
                x = x[x > 2]
                x = np.log(x)
            hist = np.log(hist)
            if(grafico == 1):
                axs[k].scatter(x,hist,c = color)
                M[i][k],A[i][k],r2 = LM(x,hist,axs[k],color)
                axs[k].grid()
            else:
                M[i][k],A[i][k],r2 = LM(x,hist,-1,color)
                
            k += 1
        if(grafico == 1):
            fig.suptitle(f"Faixa {int(i)+1}")
            fig.legend()
            plt.show()
    return -M,A

def filter_age(contacts,data,name):
    age = np.array([contacts["cnt_age_exact"].values,contacts["part_id"]]).T
    min = contacts["cnt_age_est_min"].values
    max = contacts["cnt_age_est_max"].values
    for i,j,k in zip(range(len(age)),min,max):
        try:
            age[i][0] = float(age[i][0])
            if(math.isnan(float(age[i][0]))):
                age[i][0] = (j+k)/2
        except:
            age[i][0] = -1
    age = age[age[:,0] != -1]

    ids = data["part_id"].values
    mapeamento_ids = {}
    for i in zip(ids,range(data.shape[0])):
        mapeamento_ids[i[0]] = name+"_"+str(i[1])
    
    age = pd.DataFrame(data=age, columns=["Idade","part_id"])
    age = pd.merge(age, data, on='part_id')
    age['part_id'] = age['part_id'].replace(mapeamento_ids)
    return age.values.tolist()

def get_all_data():
    dir = "./input/"
    dir = [i for i in os.listdir(dir) if("." not in i)]
    a = []
    for j in dir:

        files = os.listdir("./input/"+j+"/")
        data = ["./input/"+j+"/"+i for i in files if( 'participant' in i)][0]
        contacts = ["./input/"+j+"/"+i for i in files if( 'contact' in i)][0]
        data = pd.read_csv(data,usecols = ["part_id","part_age"])
        contacts = pd.read_csv(contacts,usecols = ["part_id","cnt_age_exact","cnt_age_est_min","cnt_age_est_max"])
        a += filter_age(contacts,data,j)
    a = pd.DataFrame(data=a, columns=["Contato_idade","id","Idade"])
    a = a.dropna(subset=['Idade'])
    return a

def transform_faixa(a,coluna,quartis):
    faixas = []
    for i in a[coluna]:
        for j in range(len(quartis)):
            if(i<quartis[j]):
                break
        faixas.append(j-1)
    faixas = np.array(faixas)
    a[coluna+"Faixas"] = faixas
    return a

def generate_distribution_byfaixas(contagem,faixas):
    graus = np.sum(contagem,axis = 1)
    A = []
    B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
    B = [np.mean(B[faixas == i,:],axis = 0) for i in range(5)]
    C = []
    S = []
    eixos = ([0,0],[0,1],[1,0],[1,1],[2,0])
    fig, axs = plt.subplots(3, 2, figsize=(10, 15))
    for faixa,eix in zip(range(5),eixos):
        g = graus[faixas == faixa]
        x = np.arange(59+1)
        y = np.zeros(len(x))
        for i in g:
            try:
                y[int(i)] += 1
            except:
                pass
        id0 = y == 0
        y = y/np.sum(y)
        soma = np.cumsum(y)
        soma[id0] = -1
        S.append(np.cumsum(y))
        #np.savetxt(f"./C/dados/distribution_{faixa}.txt",soma,fmt = "%f");
        y =y[x>2]
        x = x[x>2]
        x = x[y>0]
        y =y[y>0]
        r = 0
        max_ = 0
        for max in range(100,30,-1):
            teste_y =y[x<max]
            teste_x = x[x<max]
            angular, beta,r2 = LM(teste_x,np.log(teste_y),-1)
            if(r2 > r):
                r = r2
                max_ = max
        y =y[x<max_]
        x = x[x<max_]
        if(faixa+1 >= 4):
            x[(x < 25) & (x >15)] = np.mean(x[(x < 25) & (x >15)])
            y[(x < 25) & (x >15)] = np.mean(y[(x < 25) & (x >15)])
            x[x > 25] = np.mean(x[x>25])
            y[x > 25] = np.mean(y[x>25])
        angular, beta,r2 = LM(x,np.log(y),-1)
        C.append(beta)
        A.append(angular)

        axs[eix[0],eix[1]].scatter(x,np.log(y),alpha = 0.8,edgecolors='black', linewidths=1.)
        axs[eix[0],eix[1]].plot(x,x*angular + beta,c = 'red',label = 'R² = {:.3f}; {:.3f}x + {:.3f}'.format(r2,angular,beta))
        axs[eix[0],eix[1]].set_title(f"Faixa {faixa+1}")
        axs[eix[0],eix[1]].set_xlabel(f"Número de Ligações")
        axs[eix[0],eix[1]].set_ylabel(f"log(P)")
        axs[eix[0],eix[1]].grid()
        axs[eix[0],eix[1]].legend()
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
    plt.tight_layout()
    axs[2, 1].axis('off')
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    axs[2, 0].set_position([0.25, 0.05, 0.5, 0.25])
    plt.savefig(f"./img/contatos_faixa.png")
    plt.show()

    return -np.array(A),B,C,S

def comparacao(graus,faixas,tem_plot = 1,Graus = np.loadtxt("./C/dados/graus_finais.txt"),faixas_ = np.loadtxt("./C/dados/faixas_finais.txt")):
    M = np.zeros((5,5))
    for faixa in range(5):
        if(tem_plot == 1):
            fig, axs = plt.subplots(5,figsize=(8, 20))
        for plotinho in range(5):

            g = Graus.T[plotinho][faixas_ == faixa]
            c = graus.T[plotinho][faixas == faixa]
            C = np.sqrt(-0.5*np.log(0.05/2))*np.sqrt((len(g)+len(c))/(len(g)*len(c)))

            x1,y1 = np.unique(g,return_counts=True)
            x2,y2 = np.unique(c,return_counts=True)
            y1 = y1/np.sum(y1)
            y2 = y2/np.sum(y2)

            #isin = np.isin(x2,x1)
            #isin2 = np.isin(x1,x2)
            #D = np.max(np.abs(y1[isin2] - y2[isin]))


            if(tem_plot == 1):
                axs[plotinho].bar(x1,y1, label = 'Rede',zorder = 1,alpha = 0.75)
                axs[plotinho].bar(x2,y2, label = 'Dados',zorder = 0)
                axs[plotinho].grid()
                axs[plotinho].legend()
                axs[plotinho].set_xlim(-0.5,50)
            #
            #M[faixa][plotinho] = D < C
        if(tem_plot == 1):
            fig.suptitle(f"Faixa {faixa+1}")
            plt.show()
    print(M)