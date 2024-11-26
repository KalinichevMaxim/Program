import numpy as np
cimport numpy as np
cimport cython
from cython.parallel import prange, parallel
from libc.math cimport pi, exp

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)	
cdef double energy_calc(double[:] spin, int size, double ji, double field, double energy):
	energy = 0
	cdef int i
	for i in range(size - 1):
		energy += ji * spin[i] * spin[i+1] - field * spin[i]	
	energy += ji * spin[size - 1] * spin[0] - field * spin[size - 1]
	return energy

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)	
cdef void calc(double[:] result_mv, double[:,:] energy_mv, double[:,:] spin_mv, double[:,:] trial_spin_mv, double[:] energetic, double[:] temp_mv, double[:] rand_num_mv, long int[:] flip_arr_mv, int size_spin, int size_rand, double ji_mv, double field_mv, int k) nogil:
	cdef int i, j	
	energy_mv[0][k] = 0
	for i in range(size_spin - 1):
		energy_mv[0][k] += -ji_mv * spin_mv[i][k] * spin_mv[i+1][k] - field_mv * spin_mv[i][k]	
	energy_mv[0][k] += -ji_mv * spin_mv[size_spin - 1][k] * spin_mv[0][k] - field_mv * spin_mv[size_spin - 1][k]
	for i in range(size_rand - 1):
		for j in range(size_spin ):
			trial_spin_mv[j][k] = spin_mv[j][k]
		a = flip_arr_mv[i]
		trial_spin_mv[a][k] = - trial_spin_mv[a][k]
		energetic[k]  = 0 
		for j in range(size_spin - 1):
			energetic[k] += -ji_mv * trial_spin_mv[j][k] * trial_spin_mv[j+1][k] - field_mv * trial_spin_mv[j][k]	
		energetic[k] += -ji_mv * trial_spin_mv[size_spin - 1][k] * trial_spin_mv[0][k] - field_mv * trial_spin_mv[size_spin - 1][k]
			
		if energetic[k] <= energy_mv[i][k]:
			energy_mv[i+1][k] = energetic[k]
			for j in range(size_spin ):
				spin_mv[j][k] = trial_spin_mv[j][k]
				#result_mv[j][i+1][k] = spin_mv[j][k]
							
		if energetic[k] > energy_mv[i][k]:
			distr = exp(-(energetic[k] - energy_mv[i][k]) / temp_mv[k])
			if distr >= rand_num_mv[i]:
				energy_mv[i+1][k] = energetic[k]
				for j in range(size_spin ):
					spin_mv[j][k] = trial_spin_mv[j][k]					
					#result_mv[j][i+1][k] = spin_mv[j][k]
			if distr < rand_num_mv[i]:
				for j in range(size_spin ):					
					#result_mv[j][i+1][k] = spin_mv[j][k]	
					energy_mv[i+1][k] = energy_mv[i][k]
	for i in range(size_spin):
		result_mv[k]+= spin_mv[i][k]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)	
def metro(in_spin, in_tr_spin, rand_num, flip_arr, result,  ji, field, temp, temp1, energy):

	cdef double [:] result_mv = result
	cdef double [:,:] spin_mv = in_spin
	cdef double [:,:] trial_spin_mv = in_tr_spin
	cdef long int [:] flip_arr_mv = flip_arr
	cdef double [:] rand_num_mv = rand_num
	cdef double [:, :] energy_mv = energy	
	cdef double [:] temp_mv = temp 
	cdef double ji_mv = ji
	cdef double field_mv = field
	cdef double [:] energetic = temp1 
	cdef double distr 
	
	cdef int size_rand = np.size(rand_num)
	cdef int size_spin = np.size(spin_mv[:,0])
	cdef int size_temp = np.size(temp)
	cdef int a = 0
	cdef int i, j, k, t = 0
	
	with nogil:
		for k in prange(size_temp, num_threads = 4):
			calc(result_mv, energy_mv, spin_mv, trial_spin_mv, energetic, temp_mv, rand_num_mv, flip_arr_mv, size_spin, size_rand, ji_mv, field_mv, k)	
	return result, energy



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef double calc2d(double[:] result_mv, double[:,:] energy_mv, double[:,:,:] spin_mv, double[:,:,:] trial_spin_mv, double[:] energetic, double[:] temp_mv, double[:] rand_num_mv, long int[:,:] flip_arr_mv, int size_spin, int size_rand, double ji_mv, double field_mv, int k) nogil:
	cdef int i,j,l
	energy_mv[0][k] = 0
	for i in range(size_spin - 1):
		for l in range(size_spin - 1):
			energy_mv[0][k] +=  - ji_mv * spin_mv[i][l][k] * spin_mv[i+1][l][k] - ji_mv * spin_mv[i][l][k] * spin_mv[i][l+1][k] - field_mv * spin_mv[i][l][k]	
		energy_mv[0][k] += -ji_mv * spin_mv[i][size_spin - 1][k] * spin_mv[i+1][size_spin - 1][k]- ji_mv * spin_mv[i][size_spin - 1][k] * spin_mv[i][0][k] - field_mv * spin_mv[i][size_spin - 1][k]
		energy_mv[0][k] +=  - ji_mv * spin_mv[size_spin - 1][i][k] * spin_mv[0][i][k]- ji_mv * spin_mv[size_spin - 1][i][k] * spin_mv[size_spin - 1][i+1][k] - field_mv * spin_mv[size_spin - 1][i][k]
	energy_mv[0][k] +=  - ji_mv * spin_mv[size_spin - 1][size_spin - 1][k] * spin_mv[0][size_spin - 1][k]- ji_mv * spin_mv[size_spin - 1][size_spin - 1][k] * spin_mv[size_spin - 1][0][k] - field_mv * spin_mv[size_spin - 1][size_spin - 1][k]

	for i in range(size_rand - 1):
		for j in range(size_spin ):
			for l in range(size_spin ):
				trial_spin_mv[j][l][k] = spin_mv[j][l][k]
		a = flip_arr_mv[i,0]
		b = flip_arr_mv[i,1]
		trial_spin_mv[a][b][k] = - trial_spin_mv[a][b][k]

		energetic[k] = 0 
		for j in range(size_spin - 1):
			for l in range(size_spin - 1):
				energetic[k] +=  - ji_mv * spin_mv[j][l][k] * spin_mv[j+1][l][k]- ji_mv * spin_mv[j][l][k] * spin_mv[j][l+1][k] - field_mv * spin_mv[j][l][k]
					
			energetic[k] += -ji_mv * trial_spin_mv[j][size_spin - 1][k] * trial_spin_mv[j+1][size_spin - 1][k]- ji_mv * trial_spin_mv[j][size_spin - 1][k] * trial_spin_mv[j][0][k] - field_mv * trial_spin_mv[j][size_spin - 1][k]
			energetic[k] +=  - ji_mv * trial_spin_mv[size_spin - 1][j][k] * trial_spin_mv[0][j][k]- ji_mv * trial_spin_mv[size_spin - 1][j][k] * trial_spin_mv[size_spin - 1][j+1][k] - field_mv * trial_spin_mv[size_spin - 1][j][k]
			
		energetic[k] +=  - ji_mv * trial_spin_mv[size_spin - 1][size_spin - 1][k] * trial_spin_mv[0][size_spin - 1][k]- ji_mv * trial_spin_mv[size_spin - 1][size_spin - 1][k] * trial_spin_mv[size_spin - 1][0][k] - field_mv * trial_spin_mv[size_spin - 1][size_spin - 1][k]
			
		if energetic[k] <= energy_mv[i][k]:
			energy_mv[i+1][k] = energetic[k]

			for j in range(size_spin ):
				for l in range(size_spin ):
					spin_mv[j][l][k] = trial_spin_mv[j][l][k]			
		if energetic[k] > energy_mv[i][k]:
			distr = exp(-(energetic[k] - energy_mv[i][k]) / temp_mv[k])
			if distr >= rand_num_mv[i]:
				energy_mv[i+1][k] = energetic[k]

				for j in range(size_spin ):
					for l in range(size_spin ):
						spin_mv[j][l][k] = trial_spin_mv[j][l][k]					
			if distr < rand_num_mv[i]:
	
				energy_mv[i+1][k] = energy_mv[i][k]
	for j in range(size_spin ):
		for l in range(size_spin ):
			result_mv[k]+=spin_mv[j][l][k]	
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def metro_2d(in_spin, in_tr_spin, rand_num, flip_arr, result,  ji, field, temp, temp1, energy):

	cdef double [:] result_mv = result
	cdef double [:,:, :] spin_mv = in_spin
	cdef double [:,:,:] trial_spin_mv = in_tr_spin
	cdef long int [:, :] flip_arr_mv = flip_arr
	cdef double [:] rand_num_mv = rand_num
	cdef double [:, :] energy_mv = energy	
	cdef double [:] temp_mv = temp 
	cdef double ji_mv = ji
	cdef double field_mv = field
	cdef double [:] energetic = temp1 
	cdef double distr 
	
	cdef int size_rand = np.size(rand_num)
	cdef int size_spin = np.size(spin_mv[:,0,0])
	cdef int size_temp = np.size(temp)
	cdef int a, b
	cdef int i, j, k, l 
	
	with nogil:
		for k in prange(size_temp, num_threads = 4):
			calc2d( result_mv, energy_mv, spin_mv, trial_spin_mv, energetic, temp_mv, rand_num_mv, flip_arr_mv, size_spin, size_rand, ji_mv, field_mv, k)	
	return result, energy	
