parallel_lennard_jones_md
=========================

Prerequisites
 - Recent Python installation
 - Recent Qt installation with qmake
 - OpenMPI

Usage (see script for details)
 - Simply run the run.py Python script:
 
    python run.py

Usage on computer lab
 - You need to start mpd first:

    mpd &

Usage on Abel 
 - Compile program with 
    
    make -f make_abel

Custom geometries
 - See md_geometry.py and program/create_cylinder/create_cylinder.cpp
 - Good luck!
 - (In order to use md_statistics to calculate permeability, flow in z-direction is assumed)
