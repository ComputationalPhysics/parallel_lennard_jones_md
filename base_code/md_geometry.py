from md_statistics import *
from md_unit_converter import *
from math import pi

class MDGeometry:
	def __init__(self, md):
		self.md = md
		self.unit_converter = MDUnitConverter(md)
		self.md_statistics = MDStatistics(md, self.unit_converter)
	
	def create_cylinders(self, radius = 0.45, num_cylinders_per_dimension = 1):
		system_length = self.md_statistics.calculate_system_length()
		num_nodes = self.md.nodes_x*self.md.nodes_y*self.md.nodes_z
		max_length = max(system_length[0], system_length[1])
		print "System length: ", system_length
		
		radius = radius*max_length
		self.md.run_command("%s -O3 program/create_cylinder/create_cylinder.cpp -o create_cylinder" % self.md.compiler)
		self.md.run_command("./create_cylinder "+str(num_nodes)+" "+str(radius)+" "+str(num_cylinders_per_dimension)+" "+str(system_length[0])+" "+str(system_length[1])+" "+str(system_length[2]))
		volume = system_length[2]*num_cylinders_per_dimension**2*pi*radius**2

		self.md_statistics.update_new_volume(volume=volume)