#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Ball.py
#  
#  Copyright 2019 Kristian Farkas <farkas@farkas>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import numpy as np
from numpy.random import *

# particle.py
class Ball(object):
	
	def __init__(self, c_radius, s_radius , L, index):
		
		self.position = {"r" : (c_radius - s_radius/2)*(random_sample())**(1/2), "phi": 2*np.pi*random_sample(), "z" : L*random_sample() + s_radius}
		self.position_cartesian = { "x" : self.position["r"]*np.cos(self.position["phi"]), "y" : self.position["r"]*np.sin(self.position["phi"]), "z" : self.position["z"] }
		self.noi = index
		self.c_radius = c_radius
		self.s_radius = s_radius
		
	def cart_to_cylindrical(self):
		self.position["z"] = self.position_cartesian["z"]
		self.position["r"] = np.sqrt( self.position_cartesian["x"]**2 + self.position_cartesian["y"]**2 )
		if self.position_cartesian["x"] >= 0:
			self.position["phi"] = np.arcsin(self.position_cartesian["y"]/self.position["r"])
		else:
			self.position["phi"] = -np.arcsin(self.position_cartesian["y"]/self.position["r"]) + np.pi
	
	
	def position_change(self, unit_vectors, sphere_overlap, p, u):
		self.position_cartesian["x"] += -0.3*(unit_vectors[str(p)][str(u)]["x"]*sphere_overlap[str(p)][str(u)])
		self.position_cartesian["y"] += -0.3*(unit_vectors[str(p)][str(u)]["y"]*sphere_overlap[str(p)][str(u)])
		self.position_cartesian["z"] += -0.3*(unit_vectors[str(p)][str(u)]["z"]*sphere_overlap[str(p)][str(u)])
		
	def position_change_cylinder(self, sphere_wall_overlap, p):
		self.position_cartesian["x"] += -0.3*(self.position["r"]*np.cos(self.position["phi"])*sphere_wall_overlap[str(p)])
		self.position_cartesian["y"] += -0.3*(self.position["r"]*np.sin(self.position["phi"])*sphere_wall_overlap[str(p)])
		
	def position_change_bottom(self, sphere_bottom_overlap, p):
		self.position_cartesian["z"] += 0.3*(self.position["z"])*sphere_bottom_overlap[str(p)]
		
	def position_change_gravity(self, C):
		self.position_cartesian["z"] += -C*2*self.s_radius  
    
