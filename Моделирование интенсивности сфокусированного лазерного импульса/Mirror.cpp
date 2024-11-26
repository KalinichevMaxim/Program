#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <time.h>
#include <omp.h>
#include <string>
#include "Structures.h"
#include "Import.h"
#include "Field.h"
# define M_PI   3.14159265358
# define L_SPEED  299792458
# define STEPS  100
using namespace std;

void matrixAE(cfg_t cfg, coords in_vector, coords vector, double AE[3][3]) {
	double r_sq = sqrt(pow((in_vector.x - vector.x), 2) + pow((in_vector.y - vector.y), 2) + pow((in_vector.z - vector.z), 2));
	AE[0][0] = 2 * cfg.f * r_sq - in_vector.x * (in_vector.x - vector.x);
	AE[0][1] = -in_vector.y * (in_vector.x - vector.x);
	AE[1][0] = -in_vector.x * (in_vector.y - vector.y);
	AE[1][1] = 2 * cfg.f * r_sq - in_vector.y * (in_vector.y - vector.y);
	AE[2][0] = in_vector.x * (r_sq - (in_vector.z - vector.z));
	AE[2][1] = in_vector.y * (r_sq - (in_vector.z - vector.z));
	AE[0][2] = AE[1][2] = AE[2][2] = 0;
}

void matrixAH(cfg_t cfg, coords in_vector, coords vector, double AH[3][3]) {
	AH[0][0] = - in_vector.x * (in_vector.y - vector.y);
	AH[0][1] = 2*cfg.f*(in_vector.z - vector.z) - in_vector.y * (in_vector.y - vector.y);
	AH[1][0] = -2 * cfg.f * (in_vector.z - vector.z) + in_vector.x * (in_vector.x - vector.x);
	AH[1][1] = in_vector.y * (in_vector.x - vector.x);
	AH[2][0] = 2 * cfg.f * (in_vector.y - vector.y);
	AH[2][1] = 2 * cfg.f * (in_vector.x - vector.x);
	AH[0][2] = AH[1][2] = AH[2][2] = 0;
}
double derivative(double t, double h, double(*nv)(double, double, cfg_t), cfg_t cfg) {
	double d = (nv(t + h, 1.0, cfg) - nv(t-h, 1.0, cfg)) / (2.0 * h);
	return d;
}

coords rotation_vector(double phi, coords vector) {
	coords vector_new;

	vector_new.x = vector.x * cos(phi) - vector.z * sin(phi);
	vector_new.z = vector.x * sin(phi) + vector.z * cos(phi);
	vector_new.y = vector.y;

	return vector_new;
}

EM_comp_t rotation_field(double phi, EM_comp_t field) {
	EM_comp_t field_new;

	field_new.Re_Ex = field.Re_Ex * cos(phi) + field.Re_Ez * sin(phi);
	field_new.Re_Ez = -field.Re_Ex * sin(phi) + field.Re_Ez * cos(phi);
	field_new.Re_Ey = field.Re_Ey;
	field_new.Re_Hx = field.Re_Hx * cos(phi) + field.Re_Hz * sin(phi);
	field_new.Re_Hz = -field.Re_Hx * sin(phi) + field.Re_Hz * cos(phi);
	field_new.Re_Hy = field.Re_Hy;
	field_new.Im_Ex = field.Im_Ex * cos(phi) + field.Im_Ez * sin(phi);
	field_new.Im_Ez = -field.Im_Ex * sin(phi) + field.Im_Ez * cos(phi);
	field_new.Im_Ey = field.Im_Ey;
	field_new.Im_Hx = field.Im_Hx * cos(phi) + field.Im_Hz * sin(phi);
	field_new.Im_Hz = -field.Im_Hx * sin(phi) + field.Im_Hz * cos(phi);
	field_new.Im_Hy = field.Im_Hy;
	return field_new;
}

void create(cfg_t cfg, double& x, double& y) {
	x = (cfg.h + cfg.rad*(rand() * rand() % 100000 / 50000.0 - 1.0));
	y = ((rand() * rand() % 100000 / 50000.0 - 1.0)*cfg.rad);

}

EM_comp_t sum(double start, double step, int first, int end, EM_comp_t(*funct)(double, void*), void* param) {
	EM_comp_t result, res1;
	for (int i = 0; i < 12; i++) {
		*(&(result.Re_Ex) + i) = 0;
	}

	double x = start;
	int i;

	for (i = first; i < (end + 1); i++) {
		res1 = funct(x, param);
		for (int i = 0; i < 12; i++) {
			*(&(result.Re_Ex) + i) += *(&(res1.Re_Ex) + i);
		}
		x += step;
	}

	return result;

}

EM_comp_t simpson(double a, double b, EM_comp_t(*funct)(double, void*), void* param, int steps) {
	EM_comp_t result, res1;
	double dh = (b - a) / 2.0 / steps;

	result = funct(a, param);
	res1 = funct(b, param);
	for (int i = 0; i < 12; i++) {
		*(&(result.Re_Ex) + i) += *(&(res1.Re_Ex) + i);
	}

	res1 = sum(a + 2 * dh, 2 * dh, 1, steps - 1, funct, param);
	for (int i = 0; i < 12; i++) {
		*(&(result.Re_Ex) + i) += *(&(res1.Re_Ex) + i)*2;
	}

	res1 = sum(a + dh, 2 * dh, 1, steps, funct, param);
	for (int i = 0; i < 12; i++) {
		*(&(result.Re_Ex) + i) += *(&(res1.Re_Ex) + i) * 4;
		*(&(result.Re_Ex) + i) *= dh / 3.0;
	}

	return result;
}

struct params_set_int {
	double t;
	double R;
	cfg_t cfg;
	coords vector;
	EM_t(*in_field)(cfg_t, coords);
	EM_comp_t(*field)(double, cfg_t, coords, coords, EM_t(*)(cfg_t, coords));
	EM_comp_t(*make_envelope)(double, EM_comp_t, cfg_t);

};

EM_comp_t Int(double psi, void* p) {
	struct params_set_int* params = (struct params_set_int*)p;
	double t = (params->t);
	cfg_t cfg = (params->cfg);
	double R = (params->R);
	coords vector = (params->vector);

	coords in_vector;
	in_vector.x = R * cos(psi) + cfg.h;
	in_vector.y = R * sin(psi);
	double s = (in_vector.x * in_vector.x + in_vector.y * in_vector.y) * 0.25 / (cfg.f * cfg.f);
	in_vector.z = cfg.f * (s - 1);
	double r_sq = sqrt((in_vector.x - vector.x) * (in_vector.x - vector.x) + (in_vector.y - vector.y) * (in_vector.y - vector.y) + (in_vector.z - vector.z) * (in_vector.z - vector.z));
	double tau = t + (in_vector.z - r_sq + 2 * cfg.f);

	EM_comp_t func1 = params->field(t, cfg, in_vector, vector, params->in_field), func_r;

	//for (int i = 0; i < 12; i++) {
	//	if (i < 6) {
	//		*(&(func_r.Re_Ex) + i) = *(&(func1.Re_Ex) + i + 6) * log(16) * tau / (cfg.tau_FWHM * cfg.tau_FWHM);
	//	}
	//	else {
	//		*(&(func_r.Re_Ex) + i) = *(&(func1.Re_Ex) + i - 6) * log(16) * tau / (cfg.tau_FWHM * cfg.tau_FWHM);
	//	}
	//}
	EM_comp_t func;// , func_no_der = params->make_envelope(tau, func1, cfg), func_der = params->make_envelope(tau, func_r, cfg);
	//
	//for (int i = 0; i < 12; i++) {
	//	*(&(func.Re_Ex) + i) = *(&(func_no_der.Re_Ex) + i) + *(&(func_der.Re_Ex) + i);
	//}
	//cout << func.Re_Ex << " " << func_der.Re_Ex << endl;23
	func = params->make_envelope(tau, func1, cfg);
	for (int i = 0; i < 12; i++) {
		*(&(func.Re_Ex) + i) *= R;
	}
	return func;
}

EM_comp_t Ext(double R, void* p) {
	struct params_set_int* params = (struct params_set_int*)p;
	(params->R) = R;
	EM_comp_t result;

	result = simpson(0.0, 2. * M_PI, &Int, (void*)params, STEPS);

	return result;
}


EM_comp_t calculation_by_another_method(double t, 
	cfg_t cfg, coords vector, EM_t(*in_field)(cfg_t, coords), 
	EM_comp_t(*make_envelope)(double, EM_comp_t, cfg_t),
	EM_comp_t(*field)(double, cfg_t, coords, coords, EM_t(*)(cfg_t, coords))) {

	EM_comp_t result;

	struct params_set_int params_ext;

	params_ext.t = t;
	params_ext.cfg = cfg;
	params_ext.vector = vector;
	params_ext.in_field = in_field;
	params_ext.field = field;
	params_ext.make_envelope = make_envelope;

	result = simpson(0.0, cfg.rad, &Ext, (void*)&params_ext, STEPS);

	return result;
}

EM_comp_t calc_int(int N ,double t, cfg_t cfg,coords vector, EM_t(*in_field)(cfg_t, coords), EM_comp_t(*make_envelope)(double, EM_comp_t,cfg_t),
	EM_comp_t(*field)(double, cfg_t, coords , coords , EM_t(*)(cfg_t, coords))) {

	coords in_vector;
	double s , l=N;
	double tau, r_sq;
	EM_comp_t func, result;
	result.Re_Ex = result.Re_Ey = result.Re_Ez = result.Re_Hx = result.Re_Hy = result.Re_Hz = 0;
	result.Im_Ex = result.Im_Ey = result.Im_Ez = result.Im_Hx = result.Im_Hy = result.Im_Hz = 0;
	
	for (int i = 0; i <= N; i++) {
		create(cfg, in_vector.x, in_vector.y);
		s = (in_vector.x * in_vector.x + in_vector.y * in_vector.y) * 0.25 / (cfg.f * cfg.f);
		in_vector.z = cfg.f * (s - 1) ;
		r_sq = sqrt((in_vector.x - vector.x) * (in_vector.x - vector.x) + (in_vector.y - vector.y) * (in_vector.y - vector.y) + (in_vector.z - vector.z) * (in_vector.z - vector.z)) ;
		tau = t + (in_vector.z - r_sq  + 2 * cfg.f) ;
		if (pow((in_vector.x - cfg.h), 2) + in_vector.y * in_vector.y <= cfg.rad * cfg.rad) {
			func = field(t, cfg, in_vector, vector, in_field);
			//cout <<t<<" "<< in_vector.z - r_sq + 2 * cfg.f << endl;
			result = make_envelope(tau, func, cfg);
			//result.Re_Ex += make_envelope(tau, func.Re_Ex, cfg);// +func.Im_Ex * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h * tau, make_envelope, cfg);
			//result.Re_Ey += make_envelope(tau, func.Re_Ey, cfg);// +func.Im_Ey * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Re_Ez += make_envelope(tau, func.Re_Ez, cfg);// +func.Im_Ez * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Re_Hx += make_envelope(tau, func.Re_Hx, cfg);// +func.Im_Hx * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Re_Hx += make_envelope(tau, func.Re_Hy, cfg);// +func.Im_Hy * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Re_Hz += make_envelope(tau, func.Re_Hz, cfg);// +func.Im_Hz * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);

			////cout<< /*func.Re_Ex */ derivative(tau, h*tau, make_envelope, cfg) << " " << ( - 2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg)) << endl;
			//result.Im_Ex += make_envelope(tau, func.Im_Ex, cfg);// +func.Re_Ex * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Im_Ey += make_envelope(tau, func.Im_Ey, cfg);// +func.Re_Ey * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Im_Ez += make_envelope(tau, func.Im_Ez, cfg);// +func.Re_Ez * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Im_Hx += make_envelope(tau, func.Im_Hx, cfg);// +func.Re_Hx * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Im_Hy += make_envelope(tau, func.Im_Hy, cfg);// +func.Re_Hy * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
			//result.Im_Hz += make_envelope(tau, func.Im_Hz, cfg);// +func.Re_Hz * (-2.0 * (-2.0 * log(0.5)) * (tau / pow(cfg.tau_FWHM, 2)) * make_envelope(tau, 1, cfg));//derivative(tau, h*tau, make_envelope, cfg);
		}
	}
	for (int i = 0; i < 12; i++) {
		*(&(result.Re_Ex) + i) *= 0.25 / (cfg.f * M_PI) * cfg.rad * cfg.rad / l;
	}

	return result;
}
EM_comp_t Reflected_field_der(double t, cfg_t cfg, coords in_vector, coords vector, EM_t(*in_field)(cfg_t, coords)) {
	EM_comp_t field;
	EM_t incident_field = in_field(cfg, in_vector);
	double AE[3][3], AH[3][3];
	double r_sq = sqrt((in_vector.x - vector.x) * (in_vector.x - vector.x) + (in_vector.y - vector.y) * (in_vector.y - vector.y) + (in_vector.z - vector.z) * (in_vector.z - vector.z));
	double s = (in_vector.x * in_vector.x + in_vector.y * in_vector.y) * 0.25 / (cfg.f * cfg.f);
	double l = t + (cfg.f * (s - 1) - r_sq + 2 * cfg.f);
	in_vector.z = cfg.f * (s - 1);
	matrixAE(cfg, in_vector, vector, AE);
	matrixAH(cfg, in_vector, vector, AH);

	field.Re_Ex = -incident_field.Ex * AE[0][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Re_Ey = -incident_field.Ex * AE[1][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Re_Ez = -incident_field.Ex * AE[2][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Re_Hx = -incident_field.Ex * AH[0][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Re_Hy = -incident_field.Ex * AH[1][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Re_Hz = -incident_field.Ex * AH[2][0] * sin(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);

	field.Im_Ex = incident_field.Ex * AE[0][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Im_Ey = incident_field.Ex * AE[1][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Im_Ez = incident_field.Ex * AE[2][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Im_Hx = incident_field.Ex * AH[0][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Im_Hy = incident_field.Ex * AH[1][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	field.Im_Hz = incident_field.Ex * AH[2][0] * cos(l - atan(2 * l * log(4) / (cfg.tau_FWHM * cfg.tau_FWHM))) * sqrt(1 + 4 * log(4) * l * l / (pow(cfg.tau_FWHM, 4))) / (r_sq * r_sq);
	return field;
}
EM_comp_t Reflected_field(double t, cfg_t cfg, coords in_vector, coords vector, EM_t(*in_field)(cfg_t, coords)) {
	EM_comp_t field;
	EM_t incident_field = in_field(cfg, in_vector);
	double AE[3][3], AH[3][3];
	double r_sq = sqrt((in_vector.x - vector.x) * (in_vector.x - vector.x) + (in_vector.y - vector.y) * (in_vector.y - vector.y) + (in_vector.z - vector.z) * (in_vector.z - vector.z));
	double s = (in_vector.x * in_vector.x + in_vector.y * in_vector.y) * 0.25 / (cfg.f * cfg.f);
	double l = t + (cfg.f * (s - 1) - r_sq + 2 * cfg.f);
	in_vector.z = cfg.f * (s - 1);
	matrixAE(cfg, in_vector, vector, AE);
	matrixAH(cfg, in_vector, vector, AH);
	//field.Re_Ex = 1 / (4 * cfg.f* M_PI) * incident_field.Ex * (AE[0][0] * r_sq * sin(l)  + cos(l) * in_vector.x * (in_vector.x - vector.x)  ) / pow(r_sq, 3);
	//field.Re_Ey = 1 / (4 * cfg.f * M_PI) * incident_field.Ex * AE[1][0] * (sin(l) * r_sq - cos(l)) / pow(r_sq, 3);
	//field.Re_Ez = 1 / (4 * cfg.f * M_PI) * incident_field.Ex * (AE[2][0] * sin(l) * r_sq  + cos(l) * in_vector.x * (in_vector.z - vector.z)) / pow(r_sq, 3);
	//field.Re_Hx = 1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[0][0] * (sin(l) * r_sq - cos(l)) / pow(r_sq, 3);
	//field.Re_Hy = 1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[1][0] * (sin(l) * r_sq  - cos(l)) / pow(r_sq, 3);
	//field.Re_Hz = 1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[2][0] * (sin(l) * r_sq  - cos(l)) / pow(r_sq, 3);

	//field.Im_Ex = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * (AE[0][0] * r_sq * cos(l)  - sin(l) * in_vector.x * (in_vector.x - vector.x)) / pow(r_sq, 3);
	//field.Im_Ey = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * AE[1][0] * (cos(l) * r_sq  + sin(l)) / pow(r_sq, 3);
	//field.Im_Ez = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * (AE[2][0] * cos(l) * r_sq - sin(l) * in_vector.x * (in_vector.z - vector.z)) / pow(r_sq, 3);
	//field.Im_Hx = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[0][0] * (cos(l) * r_sq  + sin(l)) / pow(r_sq, 3);
	//field.Im_Hy = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[1][0] * (cos(l) * r_sq  + sin(l)) / pow(r_sq, 3);
	//field.Im_Hz = -1 / (4 * cfg.f * M_PI) * incident_field.Ex * AH[2][0] * (cos(l) * r_sq  + sin(l)) / pow(r_sq, 3);

	field.Re_Ex = -incident_field.Ex * AE[0][0] * sin(l) / (r_sq * r_sq);
	field.Re_Ey = -incident_field.Ex * AE[1][0] * sin(l) / (r_sq * r_sq);
	field.Re_Ez = -incident_field.Ex * AE[2][0] * sin(l) / (r_sq * r_sq);
	field.Re_Hx = -incident_field.Ex * AH[0][0] * sin(l) / (r_sq * r_sq);
	field.Re_Hy = -incident_field.Ex * AH[1][0] * sin(l) / (r_sq * r_sq);
	field.Re_Hz = -incident_field.Ex * AH[2][0] * sin(l) / (r_sq * r_sq);

	field.Im_Ex =  incident_field.Ex * AE[0][0] * cos(l) / (r_sq * r_sq);
	field.Im_Ey =  incident_field.Ex * AE[1][0] * cos(l) / (r_sq * r_sq);
	field.Im_Ez =  incident_field.Ex * AE[2][0] * cos(l) / (r_sq * r_sq);
	field.Im_Hx =  incident_field.Ex * AH[0][0] * cos(l) / (r_sq * r_sq);
	field.Im_Hy =  incident_field.Ex * AH[1][0] * cos(l) / (r_sq * r_sq);
	field.Im_Hz =  incident_field.Ex * AH[2][0] * cos(l) / (r_sq * r_sq);
	return field;
}

using namespace std;
int main() {
	srand(time(NULL));
	omp_set_num_threads(4);
	clock_t start, end;
	double prog_time;
	start = clock();
	cfg_t cfg;
	coords vector, vector1;
	cfg_parameters(cfg);
	EM_t in_field;
	EM_comp_t Refl_field, Refl_field_old;
	int steps = 70;
	double len, scale;
	double d = cfg.f1*cfg.f1+1.0, t = 0;
	//int N=10000;
	if (cfg.f1 == 1) {
		len=0.5;
	}
	else {
		len=1.0;
	}
	EM_t(*Func_Field)(cfg_t, coords) = initialize_field(cfg);
	EM_comp_t(*Func_Envelope)(double, EM_comp_t, cfg_t) = initialize_envelope(cfg);
	//cout << cfg.tau_FWHM;
	ofstream fout1;
	ofstream fout2;
	// PIC test
	//fout1.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_f=1_tau=7.5_angle=90_pic_im.txt)");
	//coords vector_E, vector1_E, vector_B, vector1_B;
	//EM_comp_t Refl_field_E, Refl_field_B;
	//steps = 20;
	//t = -75;
	//vector_E.x = 0;
	//vector_E.y = 0;
	//vector_E.z = 0;
	//vector_B.x = 0;
	//vector_B.y = 0;
	//vector_B.z = 0;
	//scale = 3;
	//for (int j = 0; j < 2 * steps; j++) {
	//	vector_E.x = (j - steps) * 3;// 40/15.0;
	//	vector_B.x = (j - steps) * 3 ;//40/15.0;
	//	for (int i = 0; i < 2 * steps; i++) {
	//		vector_E.y = (i - steps) * 3;//40 / 15.0;
	//		vector_B.y = (i - steps) * 3;// 40 / 15.0;
	//		for (int k = 0; k < 2 * 4 *steps; k++) {
	//			vector_E.z = t + (k - 4 * steps) * 1.8/4.0;//6.0 / 21.0;
	//			vector_B.z = t + 2 * M_PI / 80.0 + (k - 4 * steps) * 1.8/4.0;// 6.0 / 21.0;
	//			vector1_E = rotation_vector(cfg.phi, vector_E);
	//			vector1_B = rotation_vector(cfg.phi, vector_B);
	//			Refl_field_E = calculation_by_another_method(t, cfg, vector1_E, Func_Field, Func_Envelope, *Reflected_field_der);
	//			Refl_field_B = calculation_by_another_method(t + 2 * M_PI / 80.0, cfg, vector1_B, Func_Field, Func_Envelope, *Reflected_field_der);
	//			//Refl_field_E = calculation_by_another_method(t, cfg, vector_E, Func_Field, no_envelope, *Reflected_field_der);
	//			//Refl_field_B = calculation_by_another_method(t + 2 * M_PI / 80.0, cfg, vector_B, Func_Field, no_envelope, *Reflected_field_der);
	//			Refl_field_E = rotation_field(cfg.phi, Refl_field_E);
	//			Refl_field_B = rotation_field(cfg.phi, Refl_field_B);
	//			fout1 << vector_E.x << " " << vector_E.y << " " << vector_E.z << " " << Refl_field_E.Im_Ex << " " << Refl_field_E.Im_Ey << " " <<
	//				Refl_field_E.Im_Ez <<" "<< Refl_field_B.Im_Hx << " " << Refl_field_B.Im_Hy << " " << Refl_field_B.Im_Hz << endl;
	//		}
	//	}
	//}


	// Расчет ЭМ полей в плоскости XOZ  

	fout1.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_f=1_tau=7.5_angle=90_res.txt)");
	//fout2.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_ne_f=1_tau=7.5_angle=0_res.txt)");
	//for (t = -steps * (d - len); t <= steps*(d- len); t += steps*(d- len)/2.0) {
	//t = 0;// -steps * (d - len) / 2.0;
	//len = 1.25;
	for (t = -75; t <= 75; t += 75) {
	//for (t = -steps/5.0 * (d - len); t <= steps/5.0 * (d - len); t += steps * (d - len) / 10.0){
		vector.x = 0;
		vector.y = 0;
		vector.z = 0;
		for (int j = 0; j < 2 * steps; j++) {
			
			vector.x = (- 7.5 + j * 15 / (2.0 * steps)) * 2 * M_PI;//(j - steps) *7/50.0;	

			for (int k = 0; k < 2 * steps; k++) {
				//vector.z = t + (k - steps) * 7/50.0;

				if (t == -75) {
					vector.z = (-17.0 + k * 9.5 / (2.0 * steps)) * 2 * M_PI;
				}
				if (t == 0) {
					vector.z = (-4.75 + k * 9.5 / (2.0 * steps)) * 2 * M_PI;
				}
				if (t == 75) {
					vector.z = (7.5 + k * 9.5 / (2.0 * steps)) * 2 * M_PI;
				}

				vector1 = rotation_vector(cfg.phi + cfg.psi, vector);
				Refl_field = calculation_by_another_method( t, cfg, vector1, Func_Field, Func_Envelope, *Reflected_field_der);
				//Refl_field_old = calculation_by_another_method( t, cfg, vector1, Func_Field, no_envelope, *Reflected_field);
				Refl_field = rotation_field(cfg.phi + cfg.psi, Refl_field);
				//Refl_field_old = rotation_field(cfg.phi, Refl_field_old);
				fout1 << vector.x << " " << vector.z << " " << pow(Refl_field.Re_Ex, 2) + pow(Refl_field.Im_Ex, 2) << " " << pow(Refl_field.Re_Ey, 2) + pow(Refl_field.Im_Ey, 2)
					<< " " << pow(Refl_field.Re_Ez, 2) + pow(Refl_field.Im_Ez, 2) << " " << pow(Refl_field.Re_Hx, 2) + pow(Refl_field.Im_Hx, 2) << " " <<
					pow(Refl_field.Re_Hy, 2) + pow(Refl_field.Im_Hy, 2) << " " << pow(Refl_field.Re_Hz, 2) + pow(Refl_field.Im_Hz, 2) << " " << Refl_field.Re_Ex << " " << Refl_field.Re_Ey << " " << Refl_field.Re_Ez << endl;
				//fout2 << vector.x << " " << vector.z << " " << pow(Refl_field_old.Re_Ex, 2) + pow(Refl_field_old.Im_Ex, 2) << " " << pow(Refl_field_old.Re_Ey, 2) + pow(Refl_field_old.Im_Ey, 2)
				//	<< " " << pow(Refl_field_old.Re_Ez, 2) + pow(Refl_field_old.Im_Ez, 2) << " " << pow(Refl_field_old.Re_Hx, 2) + pow(Refl_field_old.Im_Hx, 2) << " " <<
				//	pow(Refl_field_old.Re_Hy, 2) + pow(Refl_field_old.Im_Hy, 2) << " " << pow(Refl_field_old.Re_Hz, 2) + pow(Refl_field_old.Im_Hz, 2) <<" " << Refl_field_old.Re_Ex << " " << Refl_field_old.Re_Ey << " " << Refl_field_old.Re_Ez << endl;
			}
		}
	}
	

	// Одномерный расчет ЭМ полей (срез)

	//steps = 100;
	//double var = 0.2;
	////fout1.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_f=1_tau=3_angle=0_XOY_2D_Gauss(1).txt)");
	////fout1.open(string("C:\\Users\\Ìàêñèì\\python\\Mirror_sp_distribution_f=5_tau=" + to_string(int(cfg.tau_FWHM)) + "_angle=0_XOY_2D_Gauss.txt").c_str());
	//fout2.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_ne_f=5_angle=0_XOY_2D_Gauss.txt)");
	////fout2.open(string("C:\\Users\\Ìàêñèì\\python\\Mirror_sp_distribution_ne_f=1_lambda=" + to_string(int(cfg.lambda * 1E9)) + "_angle=0_2D_wavelength.txt").c_str());
	////for (t = -steps * var * (d - len); t <= steps * var * (d - len); t += steps * var * (d - len) / 2.0) {
	//	vector.x = 0;
	//	vector.y = 0;
	//	vector.z = 0;
	//	for (int j = 0; j < 2 * steps; j++) {
	//		vector.x = (j - steps) * (cfg.f1*var);				 					
	//		vector1 = rotation_vector(cfg.phi, vector);
	//		Refl_field = calculation_by_another_method(t, cfg, vector1, Func_Field, Func_Envelope, *Reflected_field);
	//		Refl_field_old = calculation_by_another_method(t, cfg, vector1, Func_Field, no_envelope, *Reflected_field);
	//		//Refl_field = rotation_field(cfg.phi, Refl_field);
	//		Refl_field_old = rotation_field(cfg.phi, Refl_field_old);
	//		//fout1 << vector.x << " " << pow(Refl_field.Re_Ex, 2) + pow(Refl_field.Im_Ex, 2) << " " << pow(Refl_field.Re_Ey, 2) + pow(Refl_field.Im_Ey, 2)
	//		//		<< " " << pow(Refl_field.Re_Ez, 2) + pow(Refl_field.Im_Ez, 2) << " " << pow(Refl_field.Re_Hx, 2) + pow(Refl_field.Im_Hx, 2) << " " <<
	//		//		pow(Refl_field.Re_Hy, 2) + pow(Refl_field.Im_Hy, 2) << " " << pow(Refl_field.Re_Hz, 2) + pow(Refl_field.Im_Hz, 2) << " " << Refl_field.Re_Ex << " " << Refl_field.Re_Ey << " " << Refl_field.Re_Ez << endl;
	//		fout2 << vector.x << " " << pow(Refl_field_old.Re_Ex, 2) + pow(Refl_field_old.Im_Ex, 2) << " " << pow(Refl_field_old.Re_Ey, 2) + pow(Refl_field_old.Im_Ey, 2)
	//				<< " " << pow(Refl_field_old.Re_Ez, 2) + pow(Refl_field_old.Im_Ez, 2) << " " << pow(Refl_field_old.Re_Hx, 2) + pow(Refl_field_old.Im_Hx, 2) << " " <<
	//				pow(Refl_field_old.Re_Hy, 2) + pow(Refl_field_old.Im_Hy, 2) << " " << pow(Refl_field_old.Re_Hz, 2) + pow(Refl_field_old.Im_Hz, 2) << " " << Refl_field_old.Re_Ex << " " << Refl_field_old.Re_Ey << " " << Refl_field_old.Re_Ez << endl;
	//		
	//		
	//	}
	////}

	//Расчет полей в плоскости XOY

	//fout2.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_XOY_f=2.5_tau=30_angle=90_Gauss.txt)");
	////fout2.open(R"(C:\Users\Ìàêñèì\python\Mirror_sp_distribution_XOY_f=2.5_angle=90_Gauss.txt)");
	//for (t = -steps * (d - len); t <= steps * (d - len); t += steps * (d - len) / 2.0) {
	//	vector.x = 0;
	//	vector.y = 0;
	//	vector.z = t;  
	//	for (int j = 0; j < 2 * steps; j++) {
	//		vector.x = (j - steps) * cfg.f1; //(cfg.f1*0.25);	
	//		for (int k = 0; k < 2 * steps; k++) {
	//			vector.y = (k - steps) * cfg.f1; //(cfg.f1*0.25);
	//			vector1 = rotation_vector(cfg.phi, vector);
	//			Refl_field = calculation_by_another_method( t, cfg, vector1, Func_Field, Func_Envelope, *Reflected_field);
	//		//	Refl_field_old = calculation_by_another_method( t, cfg, vector1, Func_Field, no_envelope, *Reflected_field);
	//			Refl_field = rotation_field(cfg.phi, Refl_field);
	//		//	Refl_field_old = rotation_field(cfg.phi, Refl_field_old);
	//			fout2 << vector.x << " " << vector.y << " " << pow(Refl_field.Re_Ex, 2) + pow(Refl_field.Im_Ex, 2) << " " << pow(Refl_field.Re_Ey, 2) + pow(Refl_field.Im_Ey, 2)
	//				<< " " << pow(Refl_field.Re_Ez, 2) + pow(Refl_field.Im_Ez, 2) << " " << pow(Refl_field.Re_Hx, 2) + pow(Refl_field.Im_Hx, 2) << " " <<
	//				pow(Refl_field.Re_Hy, 2) + pow(Refl_field.Im_Hy, 2) << " " << pow(Refl_field.Re_Hz, 2) + pow(Refl_field.Im_Hz, 2) << " " << Refl_field.Re_Ex << " " << Refl_field.Re_Ey << " " << Refl_field.Re_Ez << endl;
	//		//	fout2 << vector.x << " " << vector.y << " " << pow(Refl_field_old.Re_Ex, 2) + pow(Refl_field_old.Im_Ex, 2) << " " << pow(Refl_field_old.Re_Ey, 2) + pow(Refl_field_old.Im_Ey, 2)
	//		//		<< " " << pow(Refl_field_old.Re_Ez, 2) + pow(Refl_field_old.Im_Ez, 2) << " " << pow(Refl_field_old.Re_Hx, 2) + pow(Refl_field_old.Im_Hx, 2) << " " <<
	//		//		pow(Refl_field_old.Re_Hy, 2) + pow(Refl_field_old.Im_Hy, 2) << " " << pow(Refl_field_old.Re_Hz, 2) + pow(Refl_field_old.Im_Hz, 2) <<" " << Refl_field_old.Re_Ex << " " << Refl_field_old.Re_Ey << " " << Refl_field_old.Re_Ez << endl;
	//		}
	//	}
	//}

	fout2.close();
	fout1.close();
	end = clock();
	prog_time = ((double)(end - start)) / CLOCKS_PER_SEC/60.0;
	cout <<"Programm time:"<< prog_time << " min."<< endl;
	return 0;
}
