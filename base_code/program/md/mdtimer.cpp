#include <mdtimer.h>
#include <mpi.h>
#include <system.h>

MDTimer::MDTimer() {
    t0 = MPI_Wtime();

    moving = 0;
    forces = 0;
    mpi = 0;
    system_initialize = 0;
    sampling = 0;
    thermostat = 0;
    io = 0;
}

void MDTimer::start_moving() {
    moving_t0 = MPI_Wtime();
}

void MDTimer::end_moving() {
    moving += MPI_Wtime() - moving_t0;
}

double MDTimer::fraction_moving() {
    double t1 = MPI_Wtime();
    return moving/(t1-t0);
}

void MDTimer::start_thermostat() {
    thermostat_t0 = MPI_Wtime();
}

void MDTimer::end_thermostat() {
    thermostat += MPI_Wtime() - thermostat_t0;
}

double MDTimer::fraction_thermostat() {
    double t1 = MPI_Wtime();
    return thermostat/(t1-t0);
}

void MDTimer::start_sampling() {
    sampling_t0 = MPI_Wtime();
}

void MDTimer::end_sampling() {
    sampling += MPI_Wtime() - sampling_t0;
}

double MDTimer::fraction_sampling() {
    double t1 = MPI_Wtime();
    return sampling/(t1-t0);
}

void MDTimer::start_io() {
    io_t0 = MPI_Wtime();
}

void MDTimer::end_io() {
    io += MPI_Wtime() - io_t0;
}

double MDTimer::fraction_io() {
    double t1 = MPI_Wtime();
    return io/(t1-t0);
}

void MDTimer::start_forces() {
    forces_t0 = MPI_Wtime();
}

void MDTimer::end_forces() {
    forces += MPI_Wtime() - forces_t0;
}

double MDTimer::fraction_forces() {
    double t1 = MPI_Wtime();
    return forces/(t1-t0);
}

void MDTimer::start_mpi() {
    mpi_t0 = MPI_Wtime();
}

void MDTimer::end_mpi() {
    mpi += MPI_Wtime() - mpi_t0;
}

double MDTimer::fraction_mpi() {
    double t1 = MPI_Wtime();
    return mpi/(t1-t0);
}

void MDTimer::start_system_initialize() {
    system_initialize_t0 = MPI_Wtime();
}

void MDTimer::end_system_initialize() {
    system_initialize += MPI_Wtime() - system_initialize_t0;
}

double MDTimer::fraction_system_initialize() {
    double t1 = MPI_Wtime();
    return system_initialize/(t1-t0);
}

void MDTimer::gather_all_nodes(System *system) {
    forces_global = 0;
    moving_global = 0;
    mpi_global = 0;
    system_initialize_global = 0;

    MPI_Reduce(&moving,&moving_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&forces,&forces_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&mpi,&mpi_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&system_initialize,&system_initialize_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    forces_global/=system->num_nodes;
    mpi_global/=system->num_nodes;
    moving_global/=system->num_nodes;
    system_initialize_global/=system->num_nodes;
}
