import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport pi
cdef extern from "<complex.h>" :
	double complex creal(double complex z)
	double complex conj(double complex z)
	
def cranky(wave_func, potential, test, x_grid_size, t_grid_size, al, step):
	cdef int x_grid = x_grid_size
	cdef int t_grid = t_grid_size
	cdef double complex[:,:] wave_func_mv = wave_func
	cdef double [:] potential_mv = potential
	alpha = np.zeros((x_grid  ), dtype = np.complex_)
	beta = np.zeros((x_grid  ), dtype = np.complex_)
	b = np.zeros((x_grid ), dtype = np.complex_)
	f = np.zeros((x_grid ), dtype = np.complex_)
	cdef double complex[:] alpha_mv = alpha
	cdef double complex[:] beta_mv = beta
	cdef double complex[:] B_mv = b
	cdef double complex[:] f_mv = f
	cdef double stepik = step
	cdef double eps = 1 
	cdef double complex alph = al
	
	cdef double max_err = 0
	cdef int i, j, k
	
	for i in range(t_grid - 1):
		for j in range(1, x_grid-1):
			B_mv[j] = 2.0 + alph + 2.0 * potential_mv[j]
		for j in range(1, x_grid-1):
			f_mv[j] = wave_func_mv[i][j - 1] + wave_func_mv[i][j + 1] + (alph - 2.0 - 2.0 * potential_mv[j]) * wave_func_mv[i][j]
		f_mv[1] = f_mv[1] + wave_func_mv[i + 1][0]
		f_mv[x_grid - 2] = f_mv[x_grid - 2] + wave_func_mv[i+1][x_grid-1]
		alpha_mv[1] = -1 / B_mv[1]
		for j in range(2, x_grid-1):
			alpha_mv[j] = -1 / (B_mv[j] + alpha_mv[j - 1])
		beta_mv[1] = f_mv[1] / B_mv[1]
		for j in range(2, x_grid-1):
			beta_mv[j] = (f_mv[j] + beta_mv[j - 1]) / (B_mv[j] + alpha_mv[j - 1])
		wave_func_mv[i + 1][x_grid - 1] = beta_mv[x_grid - 1] 
		for j in range(x_grid - 2, -1, -1):
			wave_func_mv[i+1][j] = beta_mv[j] - alpha_mv[j] * wave_func_mv[i+1][j+1] 

#	for i in range(0, t_grid - 1):
#		alpha_mv[0] = 1/(2/alph+2)
#		f_mv[0] = wave_func_mv[i+1][0] + wave_func_mv[i][0] + (2/alph-2)*wave_func_mv[i][1] + wave_func_mv[i][2] + 2 * potential_mv[0] * stepik**2
#		beta_mv[0] = 1/(2/alph+2)*f[0]
#		for j in range(1, x_grid - 1):
#			f_mv[j] = wave_func_mv[i][j-1] + (2/alph-2)*wave_func_mv[i][j] + wave_func_mv[i][j+1] + 2 * potential_mv[j] * stepik**2  
#		f_mv[x_grid-1] = wave_func_mv[i][x_grid - 1] + (2/alph-2)*wave_func_mv[i][x_grid - 2] + wave_func_mv[i][x_grid - 3] + wave_func_mv[i+1][x_grid - 1] + 2 * potential_mv[x_grid - 2] * stepik**2   
#		for j in range(0, x_grid - 1):
#			alpha_mv[j+1] = 1/(-alpha_mv[j] + (2/alph+2))
#			beta_mv[j+1] = 1/(-alpha_mv[j] + (2/alph+2))*(f_mv[j] + beta_mv[j]) 
#					
#		wave_func_mv[i+1][x_grid-1] = (f_mv[x_grid-1] + beta_mv[x_grid-1])/((2/alph+2) - alpha_mv[x_grid-1])
#		
#		for j in range(x_grid - 2 , -1, -1):
#			wave_func_mv[i+1][j] = alpha_mv[j+1]*wave_func_mv[i+1][j+1] + beta_mv[j+1]
		
	for i in range(0, t_grid):
		for j in range(x_grid):
			wave_func_mv[i][j] = float(creal(conj(wave_func_mv[i][j]) * wave_func_mv[i][j]))
	#print(wave_func)
	return wave_func
