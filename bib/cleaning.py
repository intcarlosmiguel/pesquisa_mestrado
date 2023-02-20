import math
import numpy as np
import pandas as pd
import json
import os


def generate_all():
    data = pd.read_excel('./input/RawData_ComesF.xlsx',"Planilha2",index_col=False)
    contatos = pd.read_excel('./input/RawData_ComesF.xlsx',"Contatos01",index_col=False)
    contatos02 = pd.read_excel('./input/RawData_ComesF.xlsx',"Contatos02",index_col=False)

    data = data.rename(columns={'enfant gardé chez assist. Mat.\nassistante maternelle ?': 'enfant gardé chez assist. Mat. assistante maternelle ?'})
    data = data.rename(columns={'enfant  gardé maison/en famille ?': 'enfant  gardé maison en famille ?'})

    contatos = contatos.rename(columns={ 'jour1mois1' : 'dia/mês'})
    contatos = contatos.rename(columns = {'Nombre de contacts saisis JOUR 1' : 'Número de Contatos'})

    contatos02 = contatos02.rename(columns={ 'jour2mois2' : 'dia/mês'})
    contatos02 = contatos02.rename(columns = {'Nombre de contacts saisis JOUR 2' : 'Número de Contatos'})

    df = change_age(data)
    df = change_mode(df)
    df = change_traches(df)
    df = change_column(df)
    df.head()

    criancas = df[df["Tipo de questionário"] == 0].reset_index()
    adultos = df[df["Tipo de questionário"] == 1].reset_index()
    adultos = remove_adultos(adultos,0)
    criancas = remove_adultos(criancas,1)

    a = [profissional(i,adultos) for i in adultos.index]
    c = [school(i,criancas) for i in criancas.index]
    writer = pd.ExcelWriter('./dados/output.xlsx')
    create_children(criancas,c,writer)
    create_adult(adultos,a,writer)
    writer.save()

    adultos["index"] = [i for i in range(adultos.index[-1]+1)]
    criancas["index"] = ['c'+str(i) for i in range(criancas.index[-1]+1)]

    writer = pd.ExcelWriter('./dados/questionario.xlsx')
    adultos.to_excel(writer, sheet_name='adultos',index= False)
    criancas.to_excel(writer, sheet_name='criancas',index= False)
    writer.save()

    idades = []
    c = 0
    a = 0
    for i in df['Q1_cQ1'].values:
        if(i>15):
            idades.append(str(a))
            a += 1
        else:
            idades.append('c' + str(c))
            c += 1

    data = {}
    data['id'] = [i for i in range(len(df['Q1_cQ1'].values))]
    data['Idade'] = df['Q1_cQ1']
    data['#Contatos01'] = contatos['Número de Contatos'].values
    data['#Contatos02'] = contatos02['Número de Contatos'].values
    data['id_idade'] = idades
    data = pd.DataFrame.from_dict(data)

    return df,contatos,contatos02,data

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

def get_contacts(k,i,contatos):
    contacts = []
    idade_media = contatos.iloc[k,i+2]
    genero = contatos.iloc[k,i+3]

    contacts.append(idade_media)
    contacts.append(genero)
    forma = list(contatos.iloc[k,i+4:i+11].values)
    
    contacts.append(forma)

    for j in range(i+11,i+14):
        contacts.append(contatos.iloc[k,j])
    return contacts

def get_json(contatos):

    number = [get_number(contatos.columns,i) for i in contatos.columns if 'age min' in i]
    contacts = {}
    for j in range(contatos.shape[0]):
        contacts[str(j)] = [get_contacts(j,i,contatos) for i in number[:contatos['Número de Contatos'][j]]]
        #contacts[str(j)] = [get_contacts(j,i,contatos) for i in number if(not math.isnan(contatos.iloc[j,i+2]))]
    with open('./dados/contatos_dia1.json', 'w') as f:
        json.dump(contacts, f)
    return contacts

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

def generate_mortalidade(df,vacina = 'nenhuma'):

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
    
    with open('./dados/mortalidade.json', 'w') as f:
        json.dump(Nmortos, f)
    return Nmortos

def where(vetor,condition):
    a = np.arange(len(vetor))
    return a[condition]

def contacts_to_df(contacts):
    cont = [[i]*len(contacts[i]) for i in contacts]
    cont = [j for i in cont for j in i]

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
