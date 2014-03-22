#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#define MAX_ATOM_NUM 10000000
char atom_type_string[][5] = {"Ar ", "H "};
#define POSITION_SCALING_FACTOR 3.405
using namespace std;

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify the number of cpus." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	double *positions = new double[3*MAX_ATOM_NUM];
	unsigned long *atom_types = new unsigned long[MAX_ATOM_NUM];
	unsigned long *atom_ids = new unsigned long[MAX_ATOM_NUM];
	
	ofstream file ("state.xyz", ios::out);
	
	ifstream **state_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"state_files/state%04d.bin",cpu);
		state_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;
	int num_atoms = 0;
	for(int cpu=0; cpu<cpus; cpu++) { 
		unsigned long num_atoms_this_cpu = 0;
		state_files[cpu]->read(reinterpret_cast<char*>(&num_atoms_this_cpu),sizeof(unsigned long));
		state_files[cpu]->read(reinterpret_cast<char*>(&positions[6*num_atoms]),6*num_atoms_this_cpu*sizeof(double));
		state_files[cpu]->read(reinterpret_cast<char*>(&atom_types[num_atoms]), num_atoms_this_cpu*sizeof(unsigned long));
		state_files[cpu]->read(reinterpret_cast<char*>(&atom_ids[num_atoms]), num_atoms_this_cpu*sizeof(unsigned long));
		num_atoms += num_atoms_this_cpu;
	}
	
	file << num_atoms << endl;
	file << "sup" << endl;
	for(int n=0; n<num_atoms; n++) {
		// Scaling factor converts from md units to Ångströms
    	file << atom_type_string[atom_types[n]] << positions[6*n+0]*POSITION_SCALING_FACTOR << " " << positions[6*n+1]*POSITION_SCALING_FACTOR << " " << positions[6*n+2]*POSITION_SCALING_FACTOR << " " << atom_ids[n] << endl;
    }

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		state_files[cpu]->close();
	}
	cout << "state.xyz written." << endl;

	return 0;
}