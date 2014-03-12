#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#define MAX_ATOM_NUM 10000000
char atom_type_string[][5] = {"Ar ", "H "};
#define POSITION_SCALING_FACTOR 3.405
using namespace std;

int main(int args, char *argv[]) {
	if(args < 3) {
		cout << "Please specify the number of cpus and timesteps." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	int timesteps = atoi(argv[2]);
	double *positions = new double[3*MAX_ATOM_NUM];
	unsigned long *atom_type = new unsigned long[MAX_ATOM_NUM];
	unsigned long *atom_ids = new unsigned long[MAX_ATOM_NUM];
	
	ofstream file ("movie.xyz", ios::out);
	
	ifstream **movie_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"movie_files/movie%04d.bin",cpu);
		movie_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;
	cout << "Will create movie with " << timesteps << " timesteps." << endl;
	for(int timestep=0;timestep<timesteps;timestep++) {
		int num_particles = 0;
		for(int cpu=0;cpu<cpus;cpu++) { 
			unsigned long N;
			movie_files[cpu]->read(reinterpret_cast<char*>(&N),sizeof(unsigned long));
			movie_files[cpu]->read(reinterpret_cast<char*>(&positions[3*num_particles]),3*N*sizeof(double));
			movie_files[cpu]->read(reinterpret_cast<char*>(&atom_type[num_particles]),N*sizeof(unsigned long));
			movie_files[cpu]->read(reinterpret_cast<char*>(&atom_ids[num_particles]),N*sizeof(unsigned long));
			num_particles += N;
		}
		
		file << num_particles << endl;
		file << "sup" << endl;
		for(int n=0;n<num_particles;n++) {
			// We return height - r(1) because system is inverted        	
        	file << atom_type_string[atom_type[n]] << positions[3*n+0]*POSITION_SCALING_FACTOR << " " << positions[3*n+1]*POSITION_SCALING_FACTOR << " " << positions[3*n+2]*POSITION_SCALING_FACTOR << " " << atom_ids[n] << endl;
        	// file << atom_type_string[atom_type[n]] << positions[3*n+0] << " " << positions[3*n+1] << " " << positions[3*n+2] << endl;
	    }

	    cout << "Wrote timestep " << timestep << endl;
	}

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		movie_files[cpu]->close();
	}
	cout << "Movie created." << endl;

	return 0;
}