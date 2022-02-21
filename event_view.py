
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


full_data = np.loadtxt("viewable.event", delimiter=',', skiprows=0)


hello_please = 1

history_numbers = full_data[:,0]
location = 0

backstep = np.diff(full_data[:,2]) - 1 # zero except at transitions

def scatter_set(start, data, change_set):
	return_set = [data[start, 3:]]
	start += 1
	while (not(change_set[start - 1])):
		return_set = np.append(return_set,[data[start, 3:]], axis = 0)
		start += 1
	return return_set, start

while True:
	data, location = scatter_set(location, full_data, backstep)


	print(history_numbers[location])
	# print(data)




	fig = plt.figure()
	# ax = fig.add_subplot(projection='3d')
	ay = fig.add_subplot(projection='3d')

	# ax.scatter3D(data[:,0], data[:,1], data[:,2], c=data[:,3], cmap='hot')

	ay.scatter3D(data[:,0], data[:,1], data[:,2], c=data[:,3], cmap='binary')
	ay.plot3D(data[:,0], data[:,1], data[:,2], 'black')

	plt.show()





#update