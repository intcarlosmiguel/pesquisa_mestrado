import numpy as np
import matplotlib.pyplot as plt
import os

a = np.loadtxt('./C/output/infect/800.txt',delimiter= ' ').T

#X = 300
x = np.cumsum(a[0])
#x = x[x<=X]
#plt.plot(x,a[1][:len(x)]/2029,label = 'Suscetível')
#plt.plot(x,a[2][:len(x)]/2029,label = 'Expostos')
#plt.plot(x,a[-2][:len(x)]/2029,label = 'Recuperados')
#plt.plot(x,a[-3]/2029,label = 'Hospitalizados')
#plt.plot(x,a[-1]/2029,label = 'Mortos')
plt.scatter(a[-2]/2029,a[-1]/2029,label = 'Mortos')

plt.legend()
plt.grid()
#plt.xlim(0,100)
plt.show()