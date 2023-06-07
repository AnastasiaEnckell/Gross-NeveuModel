# This program plots a heatmap and a 3D plot of the probability density evolution with the Jordan-Wigner transformation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import seaborn as sns
import pandas as pd

numb = 16 #number of qubits


with open("jordan-wigner.txt") as f:
	lines = f.readlines()

def write(lines): #write the data into arrays
	data = []
	x = []
	y = []
	pos = []
	posx = []
	index = []
	prob = []
	for i in range(0,len(lines)):
		s = lines[i].split('\n')
		data.append(s)
	for i in range(0, len(data)):
		x.append(data[i][0])
	for i in range(0, len(x)):
		y.append(x[i].split('\t'))
		pos.append(y[i][0].translate({ord(i): None for i in '[]'}))
		pos[i] = np.array(pos[i].split(','))
		if len(pos[i]) == 1:
			index.append(i)

		if len(pos[i])!= 1:
			prob.append(float(y[i][1]))
			posx.append(pos[i].astype(np.float))
	finaldata = []
	helping = np.zeros(numb)

	for i in range(0, len(posx)):
		posx[i] = [element * prob[i] for element in posx[i]]

	for i in range(0, len(index)-1):
		ind1 = index[i]-i
		ind2 = index[i+1]-i-1
		helping = np.zeros(numb)
		for j in range(ind1, ind2):
			helping = np.add(helping, posx[j])
		finaldata.append(helping)
	return finaldata

dataf = write(lines)

data_array = np.array(dataf)

x_data, y_data = np.meshgrid( np.arange(data_array.shape[1]), np.arange(data_array.shape[0]) )



# Heat map 

fig = plt.figure()

df = pd.DataFrame(data_array)

ax = sns.heatmap(df, cmap="Blues_r")
ax.set_xlabel('x', fontsize=13)
ax.set_ylabel('t', fontsize=13)


plt.savefig("jordan-wigner_heatmap.png")



# 3D plot


fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
x_data = x_data.flatten() # number of sites
y_data = y_data.flatten()
z_data = data_array.flatten()
x_data1 = x_data[::2]
y_data1 = y_data[::2]
x_data2 = x_data[1::2]
y_data2 = y_data[1::2]
z_data1 = z_data[::2]
z_data2 = z_data[1::2]


ax.set_xlabel('x', fontsize=13)
ax.set_ylabel('t', fontsize=13)


ax.bar3d(x_data1-0.4, y_data1, np.zeros(len(z_data1)), 1, 1, z_data1, color='magenta',alpha=0.5)
ax.bar3d(x_data2-0.4, y_data2, np.zeros(len(z_data2)), 1, 1, z_data2, color='cyan',alpha=0.5)

ax.azim = 300	
ax.elev = 30

ax.set_zlabel('ρ(x,t)', fontsize=13)

plt.savefig("jordan-wigner_3D.png")

