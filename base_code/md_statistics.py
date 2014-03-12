from mdconfig import *
from md_unit_converter import *
from math import pi

class MDStatistics():
	def __init__(self, md, unit_converter):
		self.md = md
		self.unit_converter = unit_converter

	def get_total_unit_cells(self):
		total_unit_cells_x = self.md.unit_cells_x*self.md.nodes_x
		total_unit_cells_y = self.md.unit_cells_y*self.md.nodes_y
		total_unit_cells_z = self.md.unit_cells_z*self.md.nodes_z

		return [total_unit_cells_x, total_unit_cells_y, total_unit_cells_z]

	def calculate_system_length(self):
		total_unit_cells = self.get_total_unit_cells()
		system_length_x = total_unit_cells[0]*self.md.FCC_b
		system_length_y = total_unit_cells[1]*self.md.FCC_b
		system_length_z = total_unit_cells[2]*self.md.FCC_b
		return [system_length_x, system_length_y, system_length_z]

	def get_volume(self, path = "./"):
		with open(path+'/volume.txt') as f:
			floats = map(float, f)

		volume = floats[0]
		return volume

	def get_num_free_atoms(self, path = "./"):
		with open(path+'/number_of_free_atoms.txt') as f:
			floats = map(float, f)
		num_free_atoms = floats[0]
		return num_free_atoms

	def get_density(self, path = "./"):
		volume = self.get_volume(path=path)
		num_free_atoms = self.get_num_free_atoms(path=path)
		density = num_free_atoms / volume
		return density

	def update_new_volume(self, volume, path = "./"):
		volume_file = open(path+'/volume.txt', 'w')
		volume_file.write(str(volume))
		volume_file.close()

	def get_ideal_gas_pressure(self, temperature, path="./"):
		density = self.get_density(path=path)
		pressure = density*temperature
		return pressure

	def calc_mean_free_path(self, path="./"):
		from math import sqrt
		density = self.get_density(path=path)
		diam = 1.0
		return 1.0/(sqrt(2)*pi*diam**2*density)

	def get_gravity_force(self, path="./"):
		import re
		ini_filename = path+"/md.ini"
		f = open(ini_filename,'r')
		content = f.read()
		f.close()
		regex_float = '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
		rg = re.compile('(gravity_force = )+('+regex_float+')', re.IGNORECASE)
		matches = rg.search(content)
		force = None
		if matches:
			force = float(matches.group(2))

		return force

	def calc_viscosity(self, path="./"):
		from math import sqrt, pi
		density = self.get_density(path=path)
		diam = 1.0
		return 5.0/(16.0*sigma**2)*sqrt(self.md.mass*self.md.temperature/pi)