#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
#include "structures.h"
#include "Import.h"

# define M_PI   3.14159265358
# define L_SPEED  299792458 

void cfg_parameters(cfg_t& cfg) {
    std::string text;
    std::ifstream in1(R"(C:\Users\Максим\source\repos\Mirror\Mirror\config.txt)");

    if (in1.is_open())
    {
        in1 >> text >> cfg.a0;
        in1 >> text >> cfg.phi0;
        in1 >> text >> cfg.phi;
        in1 >> text >> cfg.f;
        in1 >> text >> cfg.rad;
        in1 >> text >> cfg.tau_FWHM;
        in1 >> text >> cfg.lambda;
        in1 >> text >> cfg.Incident_field;
        in1 >> text >> cfg.Time_mode;
        in1 >> text >> cfg.psi;
        cfg.lambda *= 1E-6;
        double omega = 2 * M_PI * L_SPEED / (cfg.lambda);
        cfg.f1 = cfg.f;
        cfg.tau_FWHM = cfg.tau_FWHM* omega * 1E-15;
        cfg.phi *= M_PI / 180;
        cfg.psi *= M_PI / 180;
        cfg.rad *= (2.0 * M_PI) / (cfg.lambda);
        cfg.f = 2 * cfg.f * cfg.rad * (1 + cos(cfg.phi)) / 2.0;
        cfg.h = (2 * cfg.f * sin(cfg.phi)) / (1 + cos(cfg.phi));
    }
    in1.close();
}