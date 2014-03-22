from md_statistics import *
from md_unit_converter import *
from math import pi

class MDGeometry:
	def __init__(self, md):
		self.md = md
		self.unit_converter = MDUnitConverter(md)
		self.md_statistics = MDStatistics(md, self.unit_converter)
	
	def create_cylinder(self, radius = 0.45, remove_atoms_outside_outer_shell = False):
		system_length = self.md_statistics.calculate_system_length()
		num_nodes = self.md.nodes_x*self.md.nodes_y*self.md.nodes_z
		max_length = max(system_length[0], system_length[1])
		print "System length: ", system_length
		
		radius = radius*max_length
		outer_radius = radius+self.md.r_cut
		self.md.run_command("%s -O3 program/create_cylinder/create_cylinder.cpp -o create_cylinder" % self.md.compiler)
		self.md.run_command("./create_cylinder "+str(num_nodes)+" "+str(radius)+str(system_length[0])+" "+str(system_length[1])+" "+str(system_length[2])+" "+str(int(remove_atoms_outside_outer_shell))+" "+str(outer_radius))
		volume = system_length[2]*pi*radius**2

		self.md_statistics.update_new_volume(volume=volume)