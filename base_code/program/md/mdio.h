#pragma once
#include <fstream>

class System;
class Settings;
using namespace std;

class MDIO
{
public:
    System *system;
    Settings *settings;
    ofstream *movie_file;
    double *data;
    bool movie_file_open;
    FILE *energy_file;
    FILE *pressure_file;
    FILE *velocity_file;
    FILE *count_periodic_file;

    MDIO();
    void setup(System *system_);
    void save_state_to_file_binary();
    void save_state_to_movie_file();
    void load_state_from_file_binary();
    void finalize();
};
