import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import spincalc
from PIL import Image
from matplotlib.animation import FuncAnimation

class Metropolis:

	def __init__(self, iter_num, spin_num, temp_num, ji, field):
		self._iter_num = iter_num
		self._spin_num = spin_num
		self._temp_num = temp_num
		self._ji = ji
		self._field = field

		
	def get_heat_capacity(self, temp_arr, energy):
		heat_capacity = np.zeros(np.size(temp_arr))
		spin_num_1d = self._spin_num**2
		#spin_num_2d = self._spin_num**2
		for i in range(np.size(temp_arr)):
			#t = 0
			#for j in range(np.size(energy[:,i])):
			#	if energy[j,i] != 0:
			#		t+= 1
			heat_capacity[i] = (np.sum(energy[:,i]**2)/np.size(energy[:,i]) - np.sum(energy[:,i])**2/np.size(energy[:,i])**2)/temp_arr[i]**2 / spin_num
		return heat_capacity
		
	def get_mean_energy(self, temp_arr, energy):
		mean_energy = np.zeros(np.size(temp_arr))
		for i in range(np.size(temp_arr)):
			mean_energy[i] = np.sum(energy[:,i])/self._iter_num  
		return mean_energy
		
	def get_magnetisation(self, temp_arr, spin):
		magnetization = np.zeros(np.size(temp_arr))
		for i in range(np.size(temp_arr)):
			#magnetization[i] = np.sum(spin[:,i])/self._spin_num 
			magnetization[i] = spin[i]/self._spin_num 
		return magnetization

	def get_magnetisation_2d(self, temp_arr, spin):
		magnetization = np.zeros(np.size(temp_arr))
		for i in range(np.size(temp_arr)):
			magnetization[i] = spin[i]/self._spin_num**2  
		return magnetization
			
	def get_spin(self, in_spin, temp_arr):
		rand_num = np.random.rand(self._iter_num)
		flip = np.random.randint(self._spin_num, size = (self._iter_num), dtype = int)
		#result = np.zeros((self._spin_num, self._iter_num, self._temp_num))
		#result[:,0, :] = in_spin
		
		result = np.zeros((self._temp_num)) 
		in_spin1 = np.copy(in_spin)
		temp_arr1 = np.copy(temp_arr)
		energy = np.zeros((self._iter_num, self._temp_num))
		return spincalc.metro(in_spin, in_spin1, rand_num, flip, result, self._ji, self._field,  temp_arr, temp_arr1, energy)

	def get_spin_2d(self, in_spin, temp_arr):
		rand_num = np.random.rand(self._iter_num)
		flip = np.random.randint(self._spin_num, size = (self._iter_num, 2), dtype = int)
		result = np.zeros((self._temp_num))
		in_spin1 = np.copy(in_spin)
		temp_arr1 = np.copy(temp_arr)
		energy = np.zeros((self._iter_num, self._temp_num))
		return spincalc.metro_2d(in_spin, in_spin1, rand_num, flip, result, self._ji, self._field,  temp_arr, temp_arr1, energy)
				
	def get_plot(self, in_spin, temp_arr):
		result = self.get_spin(in_spin, temp_arr)	
		iters1 = np.arange(0, np.size(result[1][:,self._temp_num//2]), 1)
		plt.plot(iters1, result[1][:,self._temp_num//2])
		plt.show()
		
		iters = np.arange(0, self._iter_num, 1)
		spins = np.arange(0, self._spin_num, 1)
		#extent = (np.amin(iters), np.amax(iters), np.amin(spins), np.amax(spins))
		#fig, ax = plt.subplots(1)
		#fig = ax.imshow(result[0][:,:,self._temp_num//2],extent = extent, cmap='Greys')
		#ax.set_aspect(abs((extent[1] - extent[0])/(extent[3] - extent[2])))
		#plt.colorbar(fig, ax=ax) 
		#plt.show()
		
		fig1, ax = plt.subplots(3)
		mean_energy = self.get_mean_energy(temp_arr,  result[1][10 * self._spin_num:,:])
		ax[0].plot(temp_arr,mean_energy)
		ax[0].set_title('Mean energy') 
		
		heat_capacity = self.get_heat_capacity(temp_arr, result[1][10 * self._spin_num:,:])
		ax[1].plot(temp_arr, (self._ji/temp_arr)**2/(np.cosh(self._ji*temp_arr**(-1)))**2, color = "green")
		ax[1].plot(temp_arr,heat_capacity)
		ax[1].set_title('Heat Capacity')
		
		#magnetisation = self.get_magnetisation(temp_arr, result[0][:, self._iter_num - 1, :])
		magnetisation = self.get_magnetisation(temp_arr, result[0])
		ax[2].plot(temp_arr,magnetisation)
		ax[2].set_title('Magnetisation')
		plt.show() 
	
	def get_plot_2d(self, in_spin_2d, temp_arr):
		result = self.get_spin_2d(in_spin_2d, temp_arr)	
		iters1 = np.arange(0, np.size(result[1][:,self._temp_num//2]), 1)
		plt.plot(iters1, result[1][:,self._temp_num//2])
		plt.show()
		
		fig1, ax = plt.subplots(3)
		mean_energy = self.get_mean_energy(temp_arr,  result[1][0 * self._spin_num:,:])
		ax[0].plot(2*temp_arr,mean_energy)
		ax[0].set_xlim(0.,5)
		ax[0].set_title('Mean energy') 
		
		heat_capacity = self.get_heat_capacity(temp_arr, result[1][10 * self._spin_num**2:,:])
		#ax[1].plot(temp_arr, (self._ji/temp_arr)**2/(np.cosh(self._ji*temp_arr**(-1)))**2, color = "green")
		ax[1].plot(2*temp_arr,heat_capacity)
		ax[1].set_xlim(0.,5)
		ax[1].set_title('Heat Capacity')
		
		magnetisation = self.get_magnetisation_2d(temp_arr, result[0])
		ax[2].plot(2*temp_arr,magnetisation)
		ax[2].set_xlim(0.,5)
		ax[2].set_title('Magnetisation')
		plt.show()
		
temp_num = 80
spin_num = 800
iter_num = 800 * spin_num#**2
in_temp = 0.01
ji = 1
field = 0	

in_spin_cold = np.ones((spin_num, temp_num))
in_spin_hot = np.zeros((spin_num, temp_num))
temp_arr = np.linspace(in_temp, 5, temp_num)
for i in range(spin_num):
	for j in range(temp_num):
		if np.random.rand()> 0.5:
			in_spin_hot[i][j] = 1
		else:
			in_spin_hot[i][j] = -1

in_spin_cold_2d = np.ones((spin_num,spin_num, temp_num))
			
params = Metropolis(iter_num, spin_num, temp_num, ji, field)

params.get_plot(in_spin_cold, temp_arr)

#params.get_plot_2d(in_spin_cold_2d, temp_arr)

plt.show()
