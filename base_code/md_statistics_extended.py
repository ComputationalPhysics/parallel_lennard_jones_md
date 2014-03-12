from mdconfig import *
from md_statistics import *
from md_unit_converter import *
import matplotlib.pylab as pylab
from math import pi

class MDStatisticsExtended(MDStatistics):
	def __init__(self, md, unit_converter):
		MDStatistics.__init__(self, md, unit_converter)

	def get_number_flow_rate(self, path = "./"):
		filename = path+"/statistics/count_periodic.txt"
		count_periodic = pylab.loadtxt(filename)
		last_count_periodic = count_periodic[-1]
		time = self.unit_converter.time_from_si(last_count_periodic[0]*1e-15)
		count_z = last_count_periodic[3]
		number_flow_rate = count_z / time
		return number_flow_rate

	def calculate_permeability(self, path = "./"):
		system_length = self.calculate_system_length()
		area = system_length[0]*system_length[1]
		volume = self.get_volume(path=path)
		num_free_atoms = self.get_num_free_atoms(path=path)
		number_flow_rate = self.get_number_flow_rate(path=path)
		density = self.get_density(path=path)
		volume_per_atom = volume / num_free_atoms
		volumetric_flow_rate = number_flow_rate*volume_per_atom
		diam = 1.0
		viscosity = 5.0/(16.0*diam**2)*sqrt(self.md.mass*self.md.temperature/pi)
		force = self.get_gravity_force(path=path)
		permeability = volumetric_flow_rate*viscosity / (area * self.md.mass * density)
		return permeability