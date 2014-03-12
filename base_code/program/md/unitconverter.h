#pragma once
#include <math.h>

class UnitConverter
{
public:
    double m0;
    double L0;
    double L0_angstrom;
    double E0;
    double E0ev;
    double kb;
    double t0;
    double F0;
    double T0;
    double P0;
    double v0;
    double visc0;
    double diff0;

    UnitConverter();

    double pressure_to_SI(double P);
    double pressure_from_SI(double P);

    double temperature_to_SI(double T);
    double temperature_from_SI(double T);

    double mass_to_SI(double m);
    double mass_from_SI(double m);

    double length_to_SI(double L);
    double length_from_SI(double L);

    double force_to_SI(double F);
    double force_from_SI(double F);

    double energy_to_SI(double E);
    double energy_from_SI(double E);

    double energy_to_ev(double E);
    double energy_from_ev(double E);

    double time_to_SI(double t);
    double time_from_SI(double t);

    double velocity_to_SI(double v);
    double velocity_from_SI(double v);

    double viscosity_to_SI(double v);
    double viscosity_from_SI(double v);

    double diffusion_to_SI(double d);
    double diffusion_from_SI(double d);
};
