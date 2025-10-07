import numpy as np
import matplotlib.pyplot as plt
import sys

# Quando estiver na pasta post, digite: $ python openBinTKE.py "ID_SIM" 

run = '001'
file_path = f"./../CYLINDER/001/forces_{run}.dat"

data = np.loadtxt(file_path)
normalized_data_fx = data[:,1]
normalized_data_fy = data[:,2]
time = data[:,0]
plt.plot(time,normalized_data_fx,'k-')
#plt.show()
plt.plot(time,normalized_data_fy,'r-')
plt.show()
