#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Ball_cavity_numpy.py
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
	
	def __init__(self, c_radius, s_radius , L, cavity_radius, n):
		
		self.n = n
#		self.position = np.array([ s_radius + cavity_radius + (c_radius - s_radius - cavity_radius)*(random_sample())**(1/2),  2*np.pi*random_sample(),  L*random_sample() + s_radius ])
#		self.position_cartesian = ([  self.position[0]*np.cos(self.position[1]), self.position[0]*np.sin(self.position[1]), self.position[2] ])
#		self.noi = index
		self.c_radius = c_radius
		self.s_radius = s_radius
		self.cavity_radius = cavity_radius
		self.L = L
		self.random = np.array([[random_sample(), random_sample(), random_sample()] for i in range(n)])
		self.position = np.array([[(self.s_radius + self.cavity_radius + (self.c_radius - self.s_radius - self.cavity_radius)*(self.random[i][0])**(1/2))*np.cos(2*np.pi*self.random[i][1]),
		(self.s_radius + self.cavity_radius + (self.c_radius - self.s_radius - self.cavity_radius)*(self.random[i][0])**(1/2))*np.sin(2*np.pi*self.random[i][1]),
		self.L*self.random[i][2] + self.s_radius] for i in range(n)])

									

#	def position(self):
#		arr = np.empty((self.n,2,3))
#		for i in range(self.n):
#			arr[i][0][0] = self.s_radius + self.cavity_radius + (self.c_radius - self.s_radius - self.cavity_radius)*(random_sample())**(1/2)
#			arr[i][0][1] = 2*np.pi*random_sample()
#			arr[i][0][2] = self.L*random_sample() + self.s_radius
#			arr[i][1][0] = arr[i][0][0]*np.cos(arr[i][0][1]) 
#			arr[i][1][1] = arr[i][0][0]*np.sin(arr[i][0][1])
#			arr[i][1][2] = arr[i][0][2]
#			arr[i][0] = (self.s_radius + self.cavity_radius + (self.c_radius - self.s_radius - self.cavity_radius)*(random_sample())**(1/2))*np.cos(2*np.pi*random_sample())
#			arr[i][1] = (self.s_radius + self.cavity_radius + (self.c_radius - self.s_radius - self.cavity_radius)*(random_sample())**(1/2))*np.sin(2*np.pi*random_sample())
#			arr[i][2] = self.L*random_sample() + self.s_radius
#		return arr	
	

	def position_change(self, balls, sphere_overlap, unit_vectors, i, index):
#		add = -0.3*((unit_vectors(balls)[i][index])*((sphere_overlap(balls).reshape(self.n,self.n,1))[i][index])
		self.position[i] += np.sum(-0.3*((unit_vectors(balls)[i][index])*((sphere_overlap(balls).reshape(self.n,self.n,1))[i][index])), axis = 0)
	
	def position_change_bottom(self, sphere_bottom_overlap, i):
		self.position[i][2] += 0.3*(abs(self.position[i][2]))*sphere_bottom_overlap[i]
	
	def position_change_bottom_ext(self, sphere_bottom_overlap, i):
		self.position[i][2] += abs(self.position[i][2])
		
				
	def position_change_cylinder(self, sphere_wall_overlap, i):
		self.position[i][0:2] += -0.3*(self.position[i][0:2])*sphere_wall_overlap[i]
	
	def position_change_cavity(self, sphere_cavity_overlap, i):
		self.position[i][0:2] += 0.3*(self.position[i][0:2])*sphere_cavity_overlap[i]
		

	def position_change_gravity(self, C, i):
		self.position[i][2] += -C*2*self.s_radius  
    

	
	
#	def cart_to_cylindrical(self):
#		self.position[2] = self.position_cartesian[2]
#		self.position[0] = np.sqrt( self.position_cartesian[0]**2 + self.position_cartesian[1]**2 )
#		if self.position_cartesian["x"] >= 0:
#			self.position[1] = np.arcsin(self.position_cartesian[1]/self.position[0])
#		else:
#			self.position[1] = -np.arcsin(self.position_cartesian[1]/self.position[0]) + np.pi
#	
#	
#	def position_change(self, unit_vectors, sphere_overlap, p, u):
#		self.position_cartesian["x"] += -0.3*(unit_vectors[str(p)][str(u)]["x"]*sphere_overlap[str(p)][str(u)])
#		self.position_cartesian["y"] += -0.3*(unit_vectors[str(p)][str(u)]["y"]*sphere_overlap[str(p)][str(u)])
#		self.position_cartesian["z"] += -0.3*(unit_vectors[str(p)][str(u)]["z"]*sphere_overlap[str(p)][str(u)])
#		
#	def position_change_cylinder(self, sphere_wall_overlap, p):
#		self.position_cartesian["x"] += -0.3*(self.position["r"]*np.cos(self.position["phi"])*sphere_wall_overlap[str(p)])
#		self.position_cartesian["y"] += -0.3*(self.position["r"]*np.sin(self.position["phi"])*sphere_wall_overlap[str(p)])
#		
#	def position_change_bottom(self, sphere_bottom_overlap, p):
#		self.position_cartesian["z"] += 0.3*(self.position["z"])*sphere_bottom_overlap[str(p)]
#		
#	def position_change_cavity(self, sphere_cavity_overlap, p):
#		self.position_cartesian["x"] += 0.1*(self.position["r"]*np.cos(self.position["phi"])*sphere_cavity_overlap[str(p)])
#		self.position_cartesian["y"] += 0.1*(self.position["r"]*np.sin(self.position["phi"])*sphere_cavity_overlap[str(p)])
#		
#	def position_change_gravity(self, C):
#		self.position_cartesian["z"] += -C*2*self.s_radius  
    

		
