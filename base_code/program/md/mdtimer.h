#pragma once

class System;

class MDTimer
{
public:
    MDTimer();
    double t0;

    double moving_t0;
    double moving;
    double moving_global;
    void start_moving();
    void end_moving();
    double fraction_moving();

    double thermostat_t0;
    double thermostat;
    double thermostat_global;
    void start_thermostat();
    void end_thermostat();
    double fraction_thermostat();

    double sampling_t0;
    double sampling;
    double sampling_global;
    void start_sampling();
    void end_sampling();
    double fraction_sampling();

    double io_t0;
    double io;
    double io_global;
    void start_io();
    void end_io();
    double fraction_io();

    double forces_t0;
    double forces;
    double forces_global;
    void start_forces();
    void end_forces();
    double fraction_forces();

    double mpi_t0;
    double mpi;
    double mpi_global;
    void start_mpi();
    void end_mpi();
    double fraction_mpi();

    double system_initialize_t0;
    double system_initialize;
    double system_initialize_global;
    void start_system_initialize();
    void end_system_initialize();
    double fraction_system_initialize();

    void gather_all_nodes(System *system);
};
