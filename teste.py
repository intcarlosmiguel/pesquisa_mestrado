import numpy as np
import matplotlib.pyplot as plt
import os

x = np.loadtxt('./C/output/infect/infect_vacinado.txt').T
plt.scatter(x[0],x[2])
#plt.scatter(x[0],x[1])
plt.ylabel('Número de Hospitalizados')
plt.xlabel('Fração da Vacinados')
plt.show()
