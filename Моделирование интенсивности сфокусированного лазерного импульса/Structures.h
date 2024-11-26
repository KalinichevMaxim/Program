#pragma once
struct particle_t {
    double x;
    double y;
    double z;
    double px;
    double py;
    double pz;

};

struct EM_t {
    double Ex;
    double Ey;
    double Ez;
    double Hx;
    double Hy;
    double Hz;
};
struct EM_comp_t {
    double Re_Ex;
    double Re_Ey;
    double Re_Ez;
    double Re_Hx;
    double Re_Hy;
    double Re_Hz;
    double Im_Ex;
    double Im_Ey;
    double Im_Ez;
    double Im_Hx;
    double Im_Hy;
    double Im_Hz;
};

struct cfg_t {
    double f;
    double phi0;
    double a0;
    double h;
    double rad;
    double phi;
    double t;
    double tau_FWHM;
    double lambda;
    double f1;
    double psi;
    int Time_mode;
    int Incident_field;
};
struct coords {
    double x;
    double y;
    double z;
};