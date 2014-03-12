#include <thermostat.h>
#include <system.h>
#include <statisticssampler.h>
#include <settings.h>
#include <math.h>
#include <iostream>
#include <unitconverter.h>
#include <mdtimer.h>
#include <atom_types.h>

using namespace std;

Thermostat::Thermostat(double relaxation_time_)
{
    relaxation_time = relaxation_time_;
}

void Thermostat::apply(StatisticsSampler *sampler, System *system, const double &temperature, bool frozen_atoms_only) {
    system->mdtimer->start_thermostat();

    if(frozen_atoms_only) {
        double kinetic_energy = 0;
        for(int n=0; n<system->num_atoms_local;n++) {
            if(system->atom_type[n]==FROZEN) {
                double vx = system->velocities[3*n+0];
                double vy = system->velocities[3*n+1];
                double vz = system->velocities[3*n+2];
                kinetic_energy += 0.5*system->settings->mass*(vx*vx + vy*vy + vz*vz);
            }
        }

        double kinetic_energy_per_atom = kinetic_energy / system->num_atoms_frozen;
        double current_temperature = 2.0/3*kinetic_energy_per_atom;
        double berendsen_factor = sqrt(1 + system->dt/relaxation_time*(temperature/current_temperature - 1));

        for(int n=0; n<system->num_atoms_local;n++) {
            if(system->atom_type[n]==FROZEN) {
                for(short a=0;a<3;a++) system->velocities[3*n+a] *= berendsen_factor;
            }
        }
    } else {
        sampler->sample_temperature();
        double berendsen_factor = sqrt(1 + system->dt/relaxation_time*(temperature/sampler->temperature - 1));

        for(unsigned long n=0;n<system->num_atoms_local;n++) {
            for(short a=0;a<3;a++) system->velocities[3*n+a] *= berendsen_factor;
        }
    }

    system->mdtimer->end_thermostat();
}
