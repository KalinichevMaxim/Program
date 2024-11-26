import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import temporal_shred
from PIL import Image
from matplotlib.animation import FuncAnimation

class Equation_solver:

	def __init__(self, x_grid_num, t_grid_num, step_x, step_t, alph, pot_type):
		self._x_grid_num = x_grid_num
		self._t_grid_num = t_grid_num
		self._step_x = step_x
		self._step_t = step_t
		self._alpha = alph
		self._pot_type = pot_type
		
	def set_boundary_cond(self, matrix, conditions):
		
		matrix[:,0] = conditions[0]
		matrix[:,self._x_grid_num - 1] = conditions[2]
		matrix[0 ,:] = conditions[1]
		return matrix
	
		
	def cranky(self, conditions, in_temp):
		temp = np.zeros((self._t_grid_num, self._x_grid_num))+0j
		temp = self.set_boundary_cond(temp, conditions)
		a = np.zeros(self._x_grid_num) + 0j
		b = np.zeros(self._x_grid_num) + 0j
		f = np.zeros(self._x_grid_num) + 0j 
		return np.real(temporal_shred.cranky(temp,potential(self._x_grid_num, self._step_x, self._pot_type), in_temp, self._x_grid_num, self._t_grid_num, self._alpha, self._step_x))
		
	def get_plot(self, conditions, in_temp):
		
		temp = self.cranky(conditions, in_temp)	
		x = np.array([i * self._step_x for i in range(self._x_grid_num)])
		y = np.array([i* self._step_t for i in range(self._t_grid_num)])
		z = temp
		fig = go.Figure(go.Surface(x = np.array([i * self._step_x for i in range(self._x_grid_num)]), y = np.array([i* self._step_t for i in range(self._t_grid_num)]), z = temp, colorscale='Hot'))
		X, Y = np.meshgrid(x, y)
		#ax2.contour(X, Y, z, colors='black');
		fig.show()

def potential(x_grid_num, step, pot_type):
	potential = np.zeros(x_grid_num)
	if pot_type == 1:
		for i in range(x_grid_num):
			if i<x_grid_num//5:
				potential[i] =0.15
	if pot_type == 2:
		for i in range(x_grid_num):
			if i>0*x_grid_num//10 and i<x_grid_num//5:
				potential[i] =-0.5				
	if pot_type == 3:
		for i in range(x_grid_num):
			if i > x_grid_num//3 and i < x_grid_num//3 + 9 * x_grid_num//60 or  i > 2*x_grid_num//3 - 9 * x_grid_num//60 and i < 2*x_grid_num//3:
				potential[i] =-10	
	return potential	
	
	
x_grid_num = 1000
t_grid_num =  x_grid_num
step = 1
step_t = 1
alph = 1j/2*step_t/step**2
sigma = x_grid_num//30
pot_type = 3

x_0 = np.zeros(t_grid_num)
y_0_gauss = np.concatenate((np.zeros(3*x_grid_num//8), [np.exp(-((i - x_grid_num//8))**2/(sigma)**2)/np.sqrt(np.sqrt(np.pi/2)*sigma) * np.exp(1.j*i*step) for i in range(x_grid_num//4)], np.zeros(3*x_grid_num//8)))
y_0 = np.zeros(x_grid_num)
x_n = np.zeros(t_grid_num)
y_n = np.zeros(x_grid_num)


in_func = np.zeros((t_grid_num, x_grid_num))

conditions = np.array([x_0, y_0_gauss, x_n, y_n], dtype = 'complex')
params = Equation_solver(x_grid_num, t_grid_num, step, step_t, alph,pot_type)
in_func = params.set_boundary_cond(in_func, conditions)
nyah = params.cranky(conditions, in_func)
params.get_plot(conditions, in_func)

x = np.array([i * step for i in range(x_grid_num)])
t = np.array([i * step_t for i in range(t_grid_num)])
square = np.zeros(t_grid_num)
for i in range(t_grid_num):
	for j in range(x_grid_num):
		square[i] += nyah[i,j]
fig1, ax1 = plt.subplots(facecolor='white')	
ax1.plot(t, square)
fig1.show()

fig, ax = plt.subplots(facecolor='white')
#ax = fig.add_axes([0,0,1,1], aspect = 40)

ax.set_xlim(0, x_grid_num)
ax.set_ylim(-0.05, 1.5/sigma)

def update(frame):
	i = frame
	#ax.set_title(i)
	ax.plot(x, -potential(x_grid_num, step, pot_type)**2, color = 'cyan')
	pic = ax.plot(x, nyah[i,:], color = 'purple')  
	return pic
animation = FuncAnimation(fig, update, interval=20, blit = True, frames = t_grid_num-1)
#animation.save('sine_wave.gif', writer='Pillow')
plt.show()
