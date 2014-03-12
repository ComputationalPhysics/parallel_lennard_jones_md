#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <thermostat.h>
#include <cinifile.h>
#include <unitconverter.h>
#include <settings.h>
#include <mpi.h>
#include <mdio.h>
#include <mdtimer.h>
#include <iomanip>
#include <thermostat.h>

using namespace std;

int main(int args, char *argv[]) {
    int numprocs = 1, myid = 0;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double t_start = MPI_Wtime();

    Settings *settings = new Settings("md.ini");
    int num_nodes = settings->nodes_x*settings->nodes_y*settings->nodes_z;
    if(numprocs != num_nodes) {
        if(myid==0) cout << "Wrong number of processors. " << endl << "Config files says " << num_nodes << ". MPI started with " << numprocs << "." << endl;
        MPI_Finalize();
        return(0);
    }

    System *system = new System();
    system->setup(myid, settings);

    ifstream to_continue("Tocontinue");
    double t = 0;
    unsigned long steps = 0;
    to_continue >> t;
    to_continue >> steps;
    to_continue.close();
    system->t = t;
    system->t0 = t;
    system->steps = steps;

    Thermostat thermostat(settings->thermostat_relaxation_time);
    StatisticsSampler *sampler = new StatisticsSampler(system);

    for(int i=0;i<settings->timesteps;i++) {
        system->sample_statistics = settings->statistics_interval && ((system->steps) % settings->statistics_interval == 0);
        system->step();
        if(system->sample_statistics) sampler->sample();
        if(settings->thermostat_enabled) thermostat.apply(sampler,system,settings->temperature, false);
        if(settings->thermostat_frozen_enabled) thermostat.apply(sampler,system,settings->temperature,true);
        if(settings->create_movie) system->mdio->save_state_to_movie_file();
	}

    system->mdio->save_state_to_file_binary();

    if(myid==0) {
        double total_time = MPI_Wtime() - t_start;
        cout.precision(2);
        cout << endl << "Program finished after " << total_time << " seconds. Time analysis:" << endl;
        cout << fixed
             << "      System initialize : " << system->mdtimer->system_initialize << " s ( " << 100*system->mdtimer->fraction_system_initialize() << "%)" <<  endl
             << "      Force calculation : " << system->mdtimer->forces << " s ( " << 100*system->mdtimer->fraction_forces() << "%)" <<  endl
             << "      Moving            : " << system->mdtimer->moving << " s ( " << 100*system->mdtimer->fraction_moving() << "%)" <<  endl
             << "      Thermostat        : " << system->mdtimer->thermostat << " s ( " << 100*system->mdtimer->fraction_thermostat() << "%)" <<  endl
             << "      Sampling          : " << system->mdtimer->sampling << " s ( " << 100*system->mdtimer->fraction_sampling() << "%)" <<  endl
             << "      Disk IO           : " << system->mdtimer->io << " s ( " << 100*system->mdtimer->fraction_io() << "%)" <<  endl
             << "      MPI communication : " << system->mdtimer->mpi << " s ( " << 100*system->mdtimer->fraction_mpi() << "%)" <<  endl;
        cout << endl << settings->timesteps / total_time << " timesteps / second. " << endl;
        cout << system->num_atoms_all_global*settings->timesteps / (1000*total_time) << "k atom-timesteps / second. " << endl;
        cout << system->num_atoms_all_global*settings->timesteps / (1000*total_time*numprocs) << "k atom-timesteps / second (per node). " << endl;
        ofstream to_continue_write("Tocontinue");
        to_continue_write << system->t << " " << system->steps;
        to_continue_write.close();
    }

    MPI_Finalize();

    return 0;
}
