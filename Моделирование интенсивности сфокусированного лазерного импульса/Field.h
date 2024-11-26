#pragma once
EM_t(*initialize_field(cfg_t))(cfg_t, coords);
EM_comp_t(*initialize_envelope(cfg_t))(double, EM_comp_t, cfg_t);
EM_t incident_field_const(cfg_t, coords);
EM_t incident_field_Gauss(cfg_t, coords);
EM_t incident_field_Gauss4(cfg_t, coords);

EM_comp_t cos2_envelope(double, EM_comp_t, cfg_t);
EM_comp_t no_envelope(double, EM_comp_t, cfg_t);
EM_comp_t Gauss_envelope(double, EM_comp_t, cfg_t);
EM_comp_t Chirped_Gauss_envelope(double, EM_comp_t, cfg_t);
