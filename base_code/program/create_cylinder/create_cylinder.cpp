#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define MAX_ATOM_NUM 10000000
#define ARGON  0
#define FROZEN 1

using namespace std;

int main(int args, char *argv[]) {
	if(args < 5) {
		cout << "Please specify the number of cpus, radius, lx, ly" << endl;
		cout << "./create_cylinder int double int double double" << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	double radius = atof(argv[2]);
	double system_length[2];
	system_length[0] = atof(argv[3]);
	system_length[1] = atof(argv[4]);
	bool remove_atoms_outside_shell = false;
	
	double outer_radius = radius;
	if(args > 5) remove_atoms_outside_shell = atoi(argv[5]);
	if(args > 6) outer_radius = atof(argv[6]);
	cout << "Will create cylinder with" << endl;
	cout << "  radius      : " << radius << endl;
	cout << "  outer radius: " << outer_radius << endl;
	cout << "  system size : (" << system_length[0] << "," << system_length[1] << ")" << endl;
	
	char *filename = new char[100];
	double *data = new double[6*MAX_ATOM_NUM];
	unsigned long *atom_type= new unsigned long[MAX_ATOM_NUM];
	unsigned long *atom_ids = new unsigned long[MAX_ATOM_NUM];
	unsigned long num_particles;
	
	double radius_squared = radius*radius;
	double outer_radius_squared = outer_radius*outer_radius;

	int num_atoms_before = 0;
	int num_free_atoms_global = 0;
	int num_frozen_atoms_global = 0;
	int num_atoms_global = 0;
	int num_removed_atoms_global = 0;
	
	double cylinder_center_displacement_x = system_length[0];
	double cylinder_center_displacement_y = system_length[1];
	
	for(int cpu=0;cpu<cpus;cpu++) { 
		sprintf(filename,"state_files/state%04d.bin",cpu);
		ifstream state_file(filename,ios::in | ios::binary);

		state_file.read(reinterpret_cast<char*>(&num_particles),sizeof(unsigned long));
		state_file.read(reinterpret_cast<char*>(data),6*num_particles*sizeof(double));
		state_file.read(reinterpret_cast<char*>(atom_type),num_particles*sizeof(unsigned long));
		state_file.read(reinterpret_cast<char*>(atom_ids),num_particles*sizeof(unsigned long));
		unsigned long num_atoms_local = 0;
		unsigned long num_free_atoms_local = 0;
		unsigned long num_frozen_atoms_local = 0;
		unsigned long num_removed_atoms_local = 0;
		for(int i=0;i<num_particles;i++) {
			num_atoms_before++;

			double x = data[6*i+0];
			double y = data[6*i+1];
			bool removed_this_atom = false;

			double cylinder_center_x = cylinder_center_displacement_x*0.5;
            double cylinder_center_y = cylinder_center_displacement_y*0.5;
            double dx = (x - cylinder_center_x);
            double dy = (y - cylinder_center_y);
            double dr2 = dx*dx + dy*dy;

            if(dr2 < radius_squared) {
                atom_type[i] = ARGON;

                num_free_atoms_local++;
            } else if(dr2 < outer_radius_squared || !remove_atoms_outside_shell) {
            	atom_type[i] = FROZEN;
            	num_frozen_atoms_local++;
            } else if(remove_atoms_outside_shell) {
            	removed_this_atom = true;
            	num_removed_atoms_local++;
            }

            if(!removed_this_atom) {
            	data[6*num_atoms_local + 0] = data[6*i+0];
            	data[6*num_atoms_local + 1] = data[6*i+1];
            	data[6*num_atoms_local + 2] = data[6*i+2];
            	data[6*num_atoms_local + 3] = data[6*i+3];
            	data[6*num_atoms_local + 4] = data[6*i+4];
            	data[6*num_atoms_local + 5] = data[6*i+5];
            	atom_type[num_atoms_local] = atom_type[i];
            	atom_ids[num_atoms_local] = atom_ids[i];
            	num_atoms_local++;
            }
		}
		state_file.close();

		num_atoms_global += num_atoms_local;
		num_frozen_atoms_global += num_frozen_atoms_local;
		num_free_atoms_global += num_free_atoms_local;
		num_removed_atoms_global += num_removed_atoms_local;

		ofstream save_state_file(filename,ios::out | ios::binary);
		save_state_file.write(reinterpret_cast<char*>(&num_atoms_local),sizeof(unsigned long));
		save_state_file.write(reinterpret_cast<char*>(data),6*num_atoms_local*sizeof(double));
		save_state_file.write(reinterpret_cast<char*>(atom_type),num_atoms_local*sizeof(unsigned long));
		save_state_file.write(reinterpret_cast<char*>(atom_ids),num_atoms_local*sizeof(unsigned long));
		save_state_file.close();
	}

	// Write the number of free atoms
	cout << "Createad cylinder with " << num_frozen_atoms_global << " frozen atoms, " << num_free_atoms_global << " free atoms and " << num_removed_atoms_global << " removed atoms." << endl;
	cout << "We had " << num_atoms_before << " atoms before this process. " << endl;
	ofstream free_atoms_file("number_of_free_atoms.txt");
	free_atoms_file << num_free_atoms_global;
	free_atoms_file.close();

	return 0;
}