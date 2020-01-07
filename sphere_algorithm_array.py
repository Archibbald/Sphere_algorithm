#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sphere_algorithm_array.py
#  
#  Copyright 2020 Kristian Farkas <farkas@farkas>
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


import random
import numpy as np
from numpy.random import *
import Ball_cavity_numpy as b

def sphere_overlap(balls):
    sphere_overlap = 2*balls.s_radius - np.sqrt(-2*np.dot(balls.position,(balls.position).T) + np.repeat(np.sum((balls.position)*(balls.position),axis=1),balls.n,axis=0).reshape((balls.n,balls.n)) + (np.repeat(np.sum((balls.position)*(balls.position),axis=1),balls.n,axis=0).reshape((balls.n,balls.n))).T)
    return sphere_overlap



def unit_vectors(balls):
    unit_vectors = np.nan_to_num((np.transpose((np.repeat((balls.position),balls.n,axis=0)).reshape((balls.n,balls.n,3)),(1,0,2)) - (np.repeat((balls.position),balls.n,axis=0)).reshape((balls.n,balls.n,3)))/((np.sqrt(-2*np.dot(balls.position,(balls.position).T) + np.repeat(np.sum((balls.position)*(balls.position),axis=1),balls.n,axis=0).reshape((balls.n,balls.n)) + (np.repeat(np.sum((balls.position)*(balls.position),axis=1),balls.n,axis=0).reshape((balls.n,balls.n))).T)).reshape(balls.n,balls.n,1)))
    return unit_vectors


def sphere_wall_overlap(balls):
    sphere_wall_overlap = balls.s_radius - balls.c_radius + np.sqrt((balls.position.T)[0]**2 + (balls.position.T)[1]**2)
    return sphere_wall_overlap



def sphere_cavity_overlap(balls):
    sphere_cavity_overlap = balls.s_radius + balls.cavity_radius - np.sqrt((balls.position.T)[0]**2 + (balls.position.T)[1]**2)
    return sphere_cavity_overlap


def sphere_bottom_overlap(balls):
    sphere_bottom_overlap = balls.s_radius  - (balls.position.T)[2]
    return sphere_bottom_overlap



def first_iteration(balls):
    s_diameter = 2*balls.s_radius
    s_radius = balls.s_radius
    sph_ov_max = (sphere_overlap(balls).flatten()).sort()
    while (np.unique(sphere_overlap(balls))[-2] > 10**(-4)*s_diameter) or np.amax(sphere_wall_overlap(balls)) > 10**(-4)*s_diameter or np.amax(sphere_cavity_overlap(balls)) > 10**(-4)*s_diameter or np.amax(sphere_bottom_overlap(balls)) > 10**(-4)*s_diameter:
        for k in range(len(balls.position)):
            
            if sphere_bottom_overlap(balls)[k] > 0:
                if s_radius > sphere_bottom_overlap(balls)[k] > 0:
                    balls.position_change_bottom( sphere_bottom_overlap(balls), k )
                if sphere_bottom_overlap(balls)[k] > s_radius:
                    balls.position_change_bottom_ext( sphere_bottom_overlap(balls), k )
            
            
            if sphere_wall_overlap(balls)[k] > 0:
                balls.position_change_cylinder(sphere_wall_overlap(balls), k)
                

            if sphere_cavity_overlap(balls)[k] > 0:
                balls.position_change_cavity( sphere_cavity_overlap(balls), k )
                
                
            index = ((sphere_overlap(balls)[:,k] >= 0.0) & (sphere_overlap(balls)[:,k] != 2*balls.s_radius)) 
            balls.position_change(balls, sphere_overlap, unit_vectors,k,index)
        
        print(np.unique(sphere_overlap(balls))[-2], np.amax(sphere_wall_overlap(balls)), np.amax(sphere_cavity_overlap(balls)), np.amax(sphere_bottom_overlap(balls)))
            


def centre_of_mass(balls):
    N = len(balls.position)
    COM = np.array([sum((balls.position.T)[0])/N, sum((balls.position.T)[1])/N, sum((balls.position.T)[2])/N, 
                   (np.sqrt((sum((balls.position.T)[0]))**2 + (sum((balls.position.T)[1]))**2 + (sum((balls.position.T)[2]))**2))/N])
    return COM

        
        

def int_stage(balls):
    s_diameter = 2*balls.s_radius
    s_radius = balls.s_radius
    Cg = 0
    

    
    COMs = np.array([])
    
    COMs = np.append(COMs, 0)
    COMs = np.append(COMs, centre_of_mass(balls)[3])
   
    while( abs((COMs[-1] - COMs[-2])/COMs[-1]) > 10**(-4)*s_diameter):
        for k in range(100):
                        for p in range(len(balls.position)):
                            if np.unique(sphere_overlap(balls)[p])[-2] > 10**(-2)*s_diameter:
                                Cg -= 1.25*10**(-7)
                                if sphere_bottom_overlap(balls)[p] < 0:
                                    balls.position_change_gravity(Cg, p)



                                if sphere_bottom_overlap(balls)[p] > 0:
                                    if s_radius > sphere_bottom_overlap(balls)[p] > 0:
                                        balls.position_change_bottom( sphere_bottom_overlap(balls), p )
                                    if sphere_bottom_overlap(balls)[p] > s_radius:
                                        balls.position_change_bottom_ext( sphere_bottom_overlap(balls), p )


                                if sphere_wall_overlap(balls)[p] > 0:
                                    balls.position_change_cylinder(sphere_wall_overlap(balls), p)


                                if sphere_cavity_overlap(balls)[p] > 0:
                                    balls.position_change_cavity( sphere_cavity_overlap(balls), p )


                                index = ((sphere_overlap(balls)[:,p] >= 0.0) & (sphere_overlap(balls)[:,p] != 2*balls.s_radius)) 
                                balls.position_change(balls, sphere_overlap, unit_vectors, p, index)



                            if np.unique(sphere_overlap(balls)[p])[-2] < 0.5*10**(-2)*s_diameter:
                                Cg += 1.25*10**(-7)
                                if sphere_bottom_overlap(balls)[p] < 0:
                                    balls.position_change_gravity(Cg, p)



                                if sphere_bottom_overlap(balls)[p] > 0:
                                    if s_radius > sphere_bottom_overlap(balls)[p] > 0:
                                        balls.position_change_bottom( sphere_bottom_overlap(balls), p )
                                    if sphere_bottom_overlap(balls)[p] > s_radius:
                                        balls.position_change_bottom_ext( sphere_bottom_overlap(balls), p )


                                if sphere_wall_overlap(balls)[p] > 0:
                                    balls.position_change_cylinder(sphere_wall_overlap(balls), p)


                                if sphere_cavity_overlap(balls)[p] > 0:
                                    balls.position_change_cavity( sphere_cavity_overlap(balls), p )


                                index = ((sphere_overlap(balls)[:,p] >= 0.0) & (sphere_overlap(balls)[:,p] != 2*balls.s_radius)) 
                                balls.position_change(balls, sphere_overlap, unit_vectors, p, index)


                             

                            #print(sm_c)
                            #print(sm_s)
                            #print(sm_b)
                 #       int_dict[p+1] = centre_of_mass(balls)
        COMs = np.append(COMs, centre_of_mass(balls)[3])
        print(Cg, COMs)
    
    return(Cg, COMs)   
    
    



def final_phase(balls, Cg):
    
    s_diameter = 2*balls.s_radius
    s_radius = balls.s_radius

 
    
    COMs = np.array([])
    
    
    
    C = Cg
    
    while (np.unique(sphere_overlap(balls))[-2] > 10**(-4)*s_diameter) or np.amax(sphere_wall_overlap(balls)) > 10**(-4)*s_diameter or np.amax(sphere_cavity_overlap(balls)) > 10**(-4)*s_diameter or np.amax(sphere_bottom_overlap(balls)) > 10**(-4)*s_diameter:
        COMs = np.append(COMs, 0)
        COMs = np.append(COMs, centre_of_mass(balls)[3])
    
        while( abs((COMs[-1] - COMs[-2])/COMs[-1]) > s_diameter*10**(-4)):
    
            for k in range(100):
                            for p in range(len(balls.position)):
                                    if sphere_bottom_overlap(balls)[p] < 0:
                                        balls.position_change_gravity(C, p)





                                    if sphere_bottom_overlap(balls)[p] > 0:
                                        if s_radius > sphere_bottom_overlap(balls)[p] > 0:
                                            balls.position_change_bottom( sphere_bottom_overlap(balls), p )
                                        if sphere_bottom_overlap(balls)[p] > s_radius:
                                            balls.position_change_bottom_ext( sphere_bottom_overlap(balls), p )


                                    if sphere_wall_overlap(balls)[p] > 0:
                                        balls.position_change_cylinder(sphere_wall_overlap(balls), p)


                                    if sphere_cavity_overlap(balls)[p] > 0:
                                        balls.position_change_cavity( sphere_cavity_overlap(balls), p )


                                    index = ((sphere_overlap(balls)[:,p] >= 0.0) & (sphere_overlap(balls)[:,p] != 2*balls.s_radius)) 
                                    balls.position_change(balls, sphere_overlap, unit_vectors, p, index)


            COMs = np.append(COMs, centre_of_mass(balls)[3])
            print(C, COMs, abs((COMs[-1] - COMs[-2])/COMs[-1]), np.unique(sphere_overlap(balls))[-2], np.amax(sphere_wall_overlap(balls)), np.amax(sphere_cavity_overlap(balls)), np.amax(sphere_bottom_overlap(balls)))

        C = C/2.0
    
    return(C, COMs)
      
      
      


def bridge_s_s(balls):
    index = ((sphere_overlap(balls) >= 0) & (sphere_overlap(balls) != 2*balls.s_radius)) 
    bridges = []
    for k in range(len(index)):
        for p in range(k,len(index)):
            if index[k][p] == True:
                bridges.append([k,p])
    arr = np.array(balls.position[bridges])
    return arr



def bridge_s_w(balls):
    arr = np.array([])
    ll = []
    for k in range(len(sphere_wall_overlap(balls))):
        if sphere_wall_overlap(balls)[k] > 0:
            
            if balls.position[k][0] >= 0:
                phi = np.arcsin((balls.position[k][1])/np.sqrt(balls.position[k][0]**2 + balls.position[k][1]**2))
                ll.append(np.array([[balls.position[k][0], balls.position[k][1], balls.position[k][2]], [balls.position[k][0] + 1.1*balls.s_radius*np.cos(phi), balls.position[k][1] + 1.1*balls.s_radius*np.sin(phi), balls.position[k][2] ] ]) )
            if balls.position[k][0] < 0:
                phi = -np.arcsin((balls.position[k][1])/np.sqrt(balls.position[k][0]**2 + balls.position[k][1]**2)) + np.pi
                ll.append(np.array([[balls.position[k][0], balls.position[k][1], balls.position[k][2]], [balls.position[k][0] + 1.1*balls.s_radius*np.cos(phi), balls.position[k][1] + 1.1*balls.s_radius*np.sin(phi), balls.position[k][2] ] ]) )
    arr = np.asarray(ll)
    return (arr)
                    
def bridge_s_c(balls):
    arr = np.array([])
    ll = []
    for k in range(len(sphere_cavity_overlap(balls))):
        if sphere_cavity_overlap(balls)[k] > 0:
            
            if balls.position[k][0] >= 0:
                phi = np.arcsin((balls.position[k][1])/np.sqrt(balls.position[k][0]**2 + balls.position[k][1]**2))
                ll.append(np.array([[balls.position[k][0], balls.position[k][1], balls.position[k][2]], [balls.position[k][0] - 1.1*balls.s_radius*np.cos(phi), balls.position[k][1] - 1.1*balls.s_radius*np.sin(phi), balls.position[k][2] ] ]) )
            if balls.position[k][0] < 0:
                phi = -np.arcsin((balls.position[k][1])/np.sqrt(balls.position[k][0]**2 + balls.position[k][1]**2)) + np.pi
                ll.append(np.array([[balls.position[k][0], balls.position[k][1], balls.position[k][2]], [balls.position[k][0] - 1.1*balls.s_radius*np.cos(phi), balls.position[k][1] - 1.1*balls.s_radius*np.sin(phi), balls.position[k][2] ] ]) )
    arr = np.asarray(ll)
    return (arr)
                    
