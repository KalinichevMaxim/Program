#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include "Structures.h"
#include "Import.h"
# define M_PI   3.14159265358

using namespace std;

EM_t incident_field_const(cfg_t cfg, coords vector) {
    EM_t field;
    field.Ex = cfg.a0;
    field.Ey = 0;
    field.Ez = 0;
    field.Hx = 0;
    field.Hy = cfg.a0;
    field.Hz = 0;

    return field;
}
EM_t incident_field_Gauss(cfg_t cfg, coords vector) {
    EM_t field;
    double C = log(10.0);
    double a = 2.15678;
    field.Ex = cfg.a0 * a * exp(-((vector.x - cfg.h) * (vector.x - cfg.h) + vector.y * vector.y) * C / (cfg.rad * cfg.rad));
    field.Ey = 0;
    field.Ez = 0;
    field.Hx = 0;
    field.Hy = cfg.a0;
    field.Hz = 0;

    return field;
}
EM_t incident_field_Gauss4(cfg_t cfg, coords vector) {
    EM_t field;
    double C = 1.382;
    double a = 1.26;
    field.Ex = cfg.a0 * a * exp(-pow(((vector.x - cfg.h) * (vector.x - cfg.h) + vector.y * vector.y),3) * C / (pow(cfg.rad,6)));
    field.Ey = 0;
    field.Ez = 0;
    field.Hx = 0;
    field.Hy = cfg.a0;
    field.Hz = 0;

    return field;
}
EM_t(*initialize_field(cfg_t cfg))(cfg_t, coords) {
    EM_t(*Func_Field)( cfg_t, coords) = NULL;

    switch (cfg.Incident_field)
    {
    case 1:
        Func_Field = &incident_field_const;
        break;
    case 2:
        Func_Field = &incident_field_Gauss;
        break;
    case 3:
        Func_Field = &incident_field_Gauss4;
        break;
    }

    return Func_Field;
}

EM_comp_t cos2_envelope(double sigma, EM_comp_t field, cfg_t cfg) {

    EM_comp_t result;
    double x, cos2;
    x = 2 * acos(pow(2, (-0.25))) * sigma / cfg.tau_FWHM;
    cos2 = (-M_PI < 2 * x) && (2 * x < M_PI) ? pow(cos(x), 2) : 0;

    for (int i = 0; i < 12; i++) {
        *(&(result.Re_Ex) + i) = *(&(field.Re_Ex) + i) *  cos2;
    }

    return result;
}
EM_comp_t no_envelope(double sigma, EM_comp_t field, cfg_t cfg) {
    return field;
}

EM_comp_t Gauss_envelope(double sigma, EM_comp_t field, cfg_t cfg) {

    EM_comp_t result;
    double x, exp2;
    x = sqrt(-2.0 * log(0.5)) * (sigma / cfg.tau_FWHM);
    exp2 = exp(-pow(x, 2));

    for (int i = 0; i < 12; i++) {
        *(&(result.Re_Ex) + i) = *(&(field.Re_Ex) + i) * exp2;
    }
    return result;
}
EM_comp_t Gauss6_envelope(double sigma, EM_comp_t field, cfg_t cfg) {

    EM_comp_t result;
    double x, exp2;
    x = sqrt(-2.0 * log(0.5)) * (sigma / cfg.tau_FWHM);
    exp2 = exp(-pow(x, 6));

    for (int i = 0; i < 12; i++) {
        *(&(result.Re_Ex) + i) = *(&(field.Re_Ex) + i) * exp2;
    }
    return result;
}
EM_comp_t Chirped_Gauss_envelope(double sigma, EM_comp_t field, cfg_t cfg) {

    EM_comp_t result;
    double x, exp2,exp21, b=0.001;
    x = sqrt(-2.0 * log(0.5)) * (sigma / cfg.tau_FWHM);
    exp2 = exp(-pow(x, 2)) * sin(b*sigma*sigma);
    exp21 = exp(-pow(x, 2)) * cos(b * sigma * sigma);
    for (int i = 0; i < 12; i++) {
        if (i < 6){
            *(&(result.Re_Ex) + i) = *(&(field.Re_Ex) + i) * exp2;// - *(&(field.Re_Ex) + 6 + i) * exp2;
        }
        if (i >= 6) {
            *(&(result.Re_Ex) + i) = *(&(field.Re_Ex) + i) * exp21;// +*(&(field.Re_Ex) + i - 6) * exp2;
        }
    }

    return result;
}

EM_comp_t(*initialize_envelope(cfg_t cfg))(double, EM_comp_t, cfg_t) {
    EM_comp_t(*Func_Envelope)(double, EM_comp_t, cfg_t) = NULL;

    switch (cfg.Time_mode)
    {
    case 0:
        Func_Envelope = &no_envelope;
        break;
    case 1:
        Func_Envelope = &cos2_envelope;
        break;
    case 2:
        Func_Envelope = &Gauss_envelope;
        break;
    case 3:
        Func_Envelope = &Chirped_Gauss_envelope;
        break;
    case 4:
        Func_Envelope = &Gauss6_envelope;
        break;
    }

    return Func_Envelope;
}



