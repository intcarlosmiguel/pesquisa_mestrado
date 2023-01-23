import networkx as nx
import pandas as pd
import numpy as np
import json
from bib.cleaning import *
from bib.rede import *

contatos01 = pd.read_excel('./dados/last.xlsx',sheet_name='Contatos_01')
contatos02 = pd.read_excel('./dados/last.xlsx',sheet_name='Contatos_02')
data = pd.read_excel('./dados/last.xlsx',sheet_name='Pessoas')

degree = (data['#Contatos01'].values + data['#Contatos02'].values)/2
degree = np.array([i for i in degree if(i > 0)])
degree = np.ceil(degree).astype(int)

visualize(degree,1)
