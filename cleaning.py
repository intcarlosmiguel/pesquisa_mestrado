import math
import numpy as np
import pandas as pd

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
    Idade do sujeito da pesquisa
    Sexo do sujeito da pesquisa
    Relação do entrevistado com a criança <15 anos Assunto da pesquisa
    Idade do entrevistado para a criança <15 anos Assunto da pesquisa
    Sexo do entrevistado para a criança <15 anos Sujeito da pesquisa
    Número de pessoas em casa
    Departamento
    Mais elevado grau
    situação profissional
    Em qual indústria você trabalha
    A criança vai à escola?
    criança cuidada em casa/com a família?
    criança mantida em assistência. Mat.\ n babá?
    nb filhos na assistente
    o assistente acolhe crianças em idade escolar?
    mantidos em um berçário?
    quantas crianças na creche
    frequência de creche
    Número de crianças em sua classe
    criança que come na cantina?
    ir para o acampamento diurno?
    ir para o acampamento de verão durante a escola?
    indo para o acampamento de verão durante as férias?
    profissão que leva a muitos contatos?
    número de contatos profissionais
    Tem ou não tem mais de 20 contactos profissionais
    Número de alunos da turma
    Estudante que come na cantina?
    dia1mês1
    Validação de CONTROLE dia1mês1 no período de onda 1 ou 2
    diasemana1
    Férias dia1mês1
    Número de contatos inseridos DIA 1
    PESSOA DE IDADE
    modo(s) de viagem semana
    modo(s) de viagem nós
    faixa(s) etária(s) de contatos no ambiente profissional
    '''
    colunas = a.split('\n')
    colunas = [i.replace('    ','') for i in colunas]
    data.columns = colunas[:-1]
    return data
def remove_adultos(df):
    df = df.drop("Lien du répondant pour l'enfant <15 ans Sujet de l'enquête", axis=1)
    df = df.drop("Age du répondant pour l'enfant <15 ans Sujet de l'enquête", axis=1)
    df = df.drop("Sexe du répondant pour l'enfant <15 ans Sujet de l'enquête", axis=1)
    df = df.drop("Type de questionnaire", axis=1)
    df = df.drop("index", axis=1)
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
    age = ['adulto' if i==1 else 'crianca' for i in df['Type de questionnaire']]
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

    contacts = {
    }

    for z,j in zip(age2,range(contatos.shape[0])):
        contacts[str(z)] = [get_contacts(j,i,contatos) for i in number if not math.isnan(contatos.iloc[0,i])]
    return contacts