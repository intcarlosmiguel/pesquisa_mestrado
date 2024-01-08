import numpy as np
import pandas as pd
import os
from datetime import datetime

def generate_casos():
    if('casos.csv' not in os.listdir('./dados/')):

        casos = pd.read_csv(
            './input/caso_full.csv',
            usecols=['city','state','date','new_confirmed']
        )
        casos = casos.dropna()
        casos = casos.drop("city", axis=1)
        casos = casos.reset_index(drop=True)
        casos['date'] = [datetime.strptime(data, '%Y-%m-%d') for data in casos['date']]
        casos2020 = casos[casos['date'].isin([i for i in casos['date'] if(i.year == 2020)])]
        casos2020.to_csv('./dados/casos_nvacinados.csv', index=False)
        casos2021 = casos[casos['date'].isin([i for i in casos['date'] if(i.year == 2021)])]

        casos2021.to_csv('./dados/casos_vacinados.csv', index=False)
    else:
        casos2020 = pd.read_csv('./dados/casos_nvacinados.csv')
        casos2021 = pd.read_csv('./dados/casos_vacinados.csv')
        casos2020['date'] = [datetime.strptime(data, '%Y-%m-%d') for data in casos2020['date']]
        casos2021['date'] = [datetime.strptime(data, '%Y-%m-%d') for data in casos2021['date']]
    return casos2020,casos2021

def generate_interna():
    if('covid.csv' not in os.listdir('./dados/')):
        covid = pd.read_csv(
            './input/INFLUD20-27-03-2023.csv',
            delimiter = ';',
            usecols = ['NU_IDADE_N','EVOLUCAO','CLASSI_FIN','DT_INTERNA','SG_UF_INTE','DT_EVOLUCA','DT_SIN_PRI']
        )
        covid = covid[covid['CLASSI_FIN'] == 5]
        covid = covid.dropna()
        covid = covid.drop("CLASSI_FIN", axis=1)
        covid = covid.reset_index(drop=True)
        covid = covid.rename({'SG_UF_INTE':'state'}, axis='columns')
        covid = covid[covid['DT_INTERNA'].isin([i for i in covid['DT_INTERNA'].values if('2020' in i)])]
        covid.to_csv('./dados/covid.csv', index=False)
        covid['DT_INTERNA'] = [datetime.strptime(data, '%d/%m/%Y') for data in covid['DT_INTERNA']]
    else:
        covid = pd.read_csv('./dados/covid.csv')
    covid['DT_INTERNA'] = [datetime.strptime(data, '%d/%m/%Y') for data in covid['DT_INTERNA']]
    covid['DT_EVOLUCA'] = [datetime.strptime(data, '%d/%m/%Y') for data in covid['DT_EVOLUCA']]
    covid['DT_SIN_PRI'] = [datetime.strptime(data, '%d/%m/%Y') for data in covid['DT_SIN_PRI']]
    return covid

def generate_vacindados():

    if('./dados/vacinados1.csv' not in os.listdir('./dados/')):
        print('Oiiiiiiiiiii')
        vacinados = pd.read_csv(
            './input/INFLUD21-27-03-2023.csv',
            delimiter = ';',
            usecols = ['NU_IDADE_N',
                       'DT_SIN_PRI',
                       'EVOLUCAO',
                       'CLASSI_FIN',
                       'DT_INTERNA',
                       'SG_UF_INTE',
                       'DT_EVOLUCA',
                       'VACINA_COV',
                        'FAB_COV_1']
        )
        vacinados = vacinados[vacinados['CLASSI_FIN'] == 5]
        #vacinados = vacinados[vacinados['VACINA_COV'] == 1]
        #vacinados = vacinados.dropna()
        vacinados = vacinados.drop("CLASSI_FIN", axis=1)
        vacinados = vacinados.reset_index(drop=True)
        vacinados = vacinados.rename({'SG_UF_INTE':'state'}, axis='columns')
        vacinados.to_csv('./dados/vacinados.csv', index=False)
    else:
        vacinados = pd.read_csv('./vacinados.csv')
    
        vacinados['DT_INTERNA'] = [datetime.strptime(data, '%d/%m/%Y') for data in vacinados['DT_INTERNA']]
        vacinados['DT_EVOLUCA'] = [datetime.strptime(data, '%d/%m/%Y') for data in vacinados['DT_EVOLUCA']]
        vacinados['DT_SIN_PRI'] = [datetime.strptime(data, '%d/%m/%Y') for data in vacinados['DT_SIN_PRI']]
    return vacinados

def calctime(df,coluna1,coluna2):
    teste = np.array([abs((i- j).days) for i,j in zip(df[coluna1],df[coluna2])])
    return np.mean(teste)

def generate_rate1(covid,coluna1):
    estados = covid['state'].unique()
    resultado = []

    for estado in estados:
        covid_estado = covid[covid['state'] ==estado]
        data_min = min(covid_estado[coluna1])
        tempo = np.array([(data - data_min).days for data in covid_estado[coluna1].values])
        tempo = tempo[tempo > 0]

        hist,valor = np.unique(tempo, return_counts= True)

        sorting = np.argsort(hist)
        hist = hist[sorting]
        valor = valor[sorting]
        resultado_estado = np.array([valor[i]/hist[i] for i in range(len(hist) - 1)])
        resultado.append(np.mean(resultado_estado))
    resultado = np.array(resultado)
    return resultado

def generate_rate2(teste,covid,casos,coluna1,coluna2,coluna3,vacina = 0.21527777777):

    estados = covid['state'].unique()
    resultado = []

    for estado in estados:

        covid_estado = covid[covid['state'] ==estado]
        
        data_min = min(covid_estado[coluna1])
        tempo = np.array([(data - data_min).days for data in covid_estado[coluna1].values])
        tempo = tempo[tempo > 0]

        tempo,estagio_i = np.unique(tempo, return_counts= True)

        sorting = np.argsort(tempo)
        tempo = tempo[sorting]
        estagio_i = estagio_i[sorting]
        
        resultado_estado = []
        #s = []
        casos2 = casos[casos['state'] ==estado]
        tempo2 = np.array([(data - data_min).days for data in casos2[coluna2].values])
        casos2 = casos2[tempo2 > 0][coluna3]
        tempo2 = tempo2[tempo2 > 0]

        ncasos = [np.sum(casos2[tempo2 == i]) for i in tempo]
        if(teste == 0):
            resultado_estado = [estagio_i[i]/(tempo[i+1] - tempo[i])/(k*vacina) for i,k in zip(range(len(tempo) - 1),ncasos) if(k != 0)]
        if(teste == 1):
            resultado_estado = [estagio_i[i]/tempo[i]/(k*vacina) for i,k in zip(range(len(tempo) - 1),ncasos) if(k != 0)]
        
        resultado_estado = np.array(resultado_estado)
        resultado.append(np.mean(resultado_estado))
    resultado = np.array(resultado)
    return 1/np.mean(resultado)

def effectiveness_hospital(covid,pfizer):
    hist = [
            [0,20],
            [20,30],
            [30,50],
            [50,70],
            [70,100000]
        ]
    label = [f'{i[0]} <= Idade < {i[1]} ' for i in hist]

    dict = {}
    for i,j in zip(hist,label):
        H_vacina = pfizer[pfizer['NU_IDADE_N'].isin(range(i[0],i[1]+1))].shape[0]
        H = covid[covid['NU_IDADE_N'].isin(range(i[0],i[1]+1))].shape[0]
        dict[j] = H_vacina/H
    return dict