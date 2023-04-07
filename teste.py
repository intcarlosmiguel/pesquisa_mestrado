import numpy as np
import matplotlib.pyplot as plt
import os

a = np.loadtxt('./C/output/infect/1.txt',delimiter= ' ').T

#X = 300
x = np.cumsum(a[0])
#x = x[x<=X]
#plt.plot(x,a[1][:len(x)]/2029,label = 'Suscetível',c = 'blue')
#plt.plot(x,a[2][:len(x)]/2029,label = 'Expostos',c = 'orange')
#
#plt.plot(x,a[-2][:len(x)]/2029,label = 'Recuperados', c = 'black')
plt.scatter(x,a[-3]/2029,label = 'Hospitalizados')
plt.scatter(x,a[-1]/2029,label = 'Mortos')
#plt.scatter(a[-2]/2029,a[-1]/2029,label = 'Mortos')

plt.legend()
plt.grid()
#plt.xlim(0,100)
plt.show()