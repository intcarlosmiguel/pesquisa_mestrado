import numpy as np

a = np.loadtxt("./output/prob_matrix.txt")
a = a/np.sum(a,axis = 1)[:, np.newaxis]
np.savetxt("./output/matrix_prob.txt",a,fmt = "%f")