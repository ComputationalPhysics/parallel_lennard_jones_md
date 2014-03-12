from mdconfig import *
from md_unit_converter import *

program = MD(dt=0.02, name="main")
md = program.compile(skip_compile=False)
unit_converter = MDUnitConverter(program)

# nodex_i is the number of processors in the i'th dimension
# Total number of processors is nodes_x*nodes_y*nodes_z
program.nodes_x = 1
program.nodes_y = 1
program.nodes_z = 1
# Number of unit cells per cpu per dimension
program.unit_cells_x = 5
program.unit_cells_y = 5
program.unit_cells_z = 5

temperature_si = 300
temperature = unit_converter.temperature_from_si(temperature_si)

# Create new system
program.reset()
program.prepare_new_system()
program.run()

print "### Applying thermostat, T="+str(temperature_si)+"K ###"
program.prepare_thermostat(temperature=temperature, timesteps=2000, run=True, save_state_path="states/01_T_300K")

# Load this state (not necessary because the last state is always in the run directory, but this is how you do it)
program.load_state(path="states/01_T_300K")
print "### Thermalizing ... ###"
program.prepare_thermalize(timesteps=10000, run=True, save_state_path="states/02_thermalized")

# Create movie XYZ-file. If you want 1 frame, just run with 1 frame
program.create_movie_files = True
movie_frames = 1000
program.prepare_thermalize(timesteps=movie_frames, run=True)
program.create_movie(frames=movie_frames)