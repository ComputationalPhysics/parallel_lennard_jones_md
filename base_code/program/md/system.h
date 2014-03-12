#pragma once

class Atom;
class ThreadControl;
class Settings;
class MDIO;
class MDTimer;
class Random;
class UnitConverter;

#include <fstream>
#include <vector>
#define EMPTY -1

using namespace std;

class System {
private:
    void initialize();
    void move();
    void mpi_move();
    void mpi_copy();
    void calculate_accelerations();
    void calculate_accelerations_many_frozen_atoms();
    void apply_gravity();
    void full_kick();
    void half_kick();

    void set_topology();
    void init_parameters();
    void create_FCC();
    void count_frozen_atoms();
    void reset();
    inline bool atom_did_change_node(double* ri, int ku);
    inline bool atom_should_be_copied(double *ri, int ku);
    inline void cell_index_from_ijk(const int &i, const int &j, const int &k, unsigned int &cell_index);
    inline void cell_index_from_vector(unsigned int *mc, unsigned int &cell_index);
public:
    Settings *settings;
    MDIO *mdio;
    Random *rnd;
    MDTimer *mdtimer;
    UnitConverter *unit_converter;
    int  *head_all_atoms;
    int  *head_free_atoms;

    bool sample_statistics;
    unsigned long steps;
    unsigned int myid;
    unsigned int node_index[3];
    unsigned int num_processors[3];
    unsigned int neighbor_nodes[6];
    short my_parity[3];
    double cell_length[3];
    double node_length[3];
    double system_length[3];
    int num_nodes;
    int max_number_of_atoms;
    int max_number_of_cells;
    unsigned long num_atoms_local;
    unsigned long num_atoms_all_global;
    unsigned long num_atoms_free;
    unsigned long num_atoms_free_global;
    unsigned long num_atoms_frozen_global;
    unsigned long num_atoms_frozen;
    unsigned long num_atoms_ghost;


    long i,j,k,n,m,a,b,c, nx, ny, nz;
    long count_periodic[3];

    double origo[3];
    double r_cut, dt, dt_half, potential_energy, t, t0, volume;
    unsigned int mc[3];  // Usually cell index vector
    unsigned int mc1[3]; // Usually cell index vector
    unsigned int num_cells_including_ghosts_yz,cell_index, cell_index_2,num_cells_including_ghosts_xyz;
    unsigned int num_cells_local[3];
    unsigned int num_cells_including_ghosts[3];
    double dr[3];
    double shift_vector[6][3];
    unsigned int **move_queue;

    double *mpi_send_buffer;
    double *mpi_receive_buffer;
    bool *atom_moved;
    unsigned long *atom_ids;
    double *positions;
    double *accelerations;
    double *velocities;

    double mass_inverse, pressure_forces;

    unsigned long *atom_type;
    int *linked_list_all_atoms;
    int *linked_list_free_atoms;
    bool *is_ghost_cell;
    double *initial_positions;
    void apply_harmonic_oscillator();

    System();
    void setup(int myid_, Settings *settings_);
    void step();
    double get_volume() {
        return system_length[0]*system_length[1]*system_length[2];
    }
};
