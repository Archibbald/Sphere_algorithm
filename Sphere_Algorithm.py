#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Sphere_Algorithm.py
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


import random
import numpy as np
from numpy.random import *
import Ball as b



def distance_balls(balls):
    distances = {}
    for b in balls:
        distances[str(b)] = {str(c): np.sqrt(balls[b].position["r"]**2 + balls[c].position["r"]**2 
                    + balls[b].position["z"]**2 + balls[c].position["z"]**2 -2*balls[b].position["z"]*balls[c].position["z"]
                    - 2*balls[b].position["r"]*balls[c].position["r"]*(np.sin(balls[b].position["phi"])*np.sin(balls[c].position["phi"]) + 
                     np.cos(balls[b].position["phi"])*np.cos(balls[c].position["phi"]))) for c in balls if b != c}
    return distances




def sphere_overlap(distance_balls, balls):
    s_radius = balls["1"].s_radius
    overlap = {}
    for b in distance_balls:
        overlap[b] = {v : -distance_balls[b][v] + 2*s_radius for v in distance_balls[b]}
    return overlap

    
    
    
def sphere_wall_overlap(balls):
    s_radius = balls["1"].s_radius
    c_radius = balls["1"].c_radius
    s_w_overlap = {}
    for b in balls:
        s_w_overlap[b] = balls[b].position["r"] + s_radius - c_radius
    return s_w_overlap


def sphere_bottom_overlap(balls):
    s_radius = balls["1"].s_radius
    s_b_overlap = {}
    for b in balls:
        s_b_overlap[b] =  -balls[b].position["z"] + s_radius
    return s_b_overlap



def unit_vectors(balls, distance_balls):
    vectors = {}
    for b in balls:
        vectors[b] = {u : {"x" : -(balls[b].position_cartesian["x"] - balls[u].position_cartesian["x"])/distance_balls[b][u],
                       "y" : -(balls[b].position_cartesian["y"] - balls[u].position_cartesian["y"])/distance_balls[b][u],
                       "z" : -(balls[b].position_cartesian["z"] - balls[u].position_cartesian["z"])/distance_balls[b][u]} for u in distance_balls[b]}
    return vectors



def sigma_max(sphere_overlap, sphere_wall_overlap, sphere_bottom_overlap):
    sigm_m = {"sphere" : { p: max(list(sphere_overlap[str(p)].values())) for p in sphere_overlap }, 
              "cylinder" : { max(list(sphere_wall_overlap.values())) },
              "bottom" :  { max(list(sphere_bottom_overlap.values())) }}
    return sigm_m
    
    
    
def ball_dict(u, s_r, c_r, L):
    balls = {str(i+1): b.Ball(c_r, s_r, L, u+1) for i in range(u)}
    return balls
    
    

def centre_of_mass(balls):
    N = len(balls)
    COM = {}
    COM["x"] = (sum([balls[str(p)].position_cartesian["x"] for p in balls]))/N
    COM["y"] = (sum([balls[str(p)].position_cartesian["y"] for p in balls]))/N
    COM["z"] = (sum([balls[str(p)].position_cartesian["z"] for p in balls]))/N
    return COM

        
def centre_of_mass_length(COM):
    l = np.sqrt( COM["x"]**2 + COM["y"]**2 + COM["z"]**2 )
    return l



def first_iteration(balls):
    
    s_diameter = 2*balls["1"].s_radius
        
    sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
    sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
    sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]

    while (max(list(sm_s.values())) > 10**(-4)*s_diameter) or max(sm_c) > 10**(-4)*s_diameter or max(sm_b) > 10**(-4)*s_diameter:
        
        for p in sphere_overlap(distance_balls(balls), balls):
            if sphere_wall_overlap(balls)[str(p)] > 0:
                balls[str(p)].position_change_cylinder( sphere_wall_overlap(balls), p )
                balls[str(p)].cart_to_cylindrical()
            
            if sphere_bottom_overlap(balls)[str(p)] > 0:
                balls[str(p)].position_change_bottom( sphere_bottom_overlap(balls), p )
                balls[str(p)].cart_to_cylindrical()
            
            
            for u in sphere_overlap(distance_balls(balls),balls)[str(p)]:
                if sphere_overlap(distance_balls(balls),balls)[str(p)][str(u)] > 0:
                    balls[str(p)].position_change(unit_vectors(balls, distance_balls(balls)), sphere_overlap(distance_balls(balls),balls),p,u)
                    balls[str(p)].cart_to_cylindrical()
        sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
        print(sm_c)
        sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
        print(sm_s)
        sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]
        print(sm_b)






def int_stage(balls):
    s_diameter = 2*balls["1"].s_radius
    Cg = 0
    
    sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
    sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
    sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]
    
    
    COMs = np.array([])
    
    COMs = np.append(COMs, 0)
    COMs = np.append(COMs, centre_of_mass_length(centre_of_mass(balls)))
   
       
    
            
#    while max(list(sm_s.values())) < 0.5*10**(-2) or max(list(sm_s.values())) > 10**(-2):
#    for p in range(10)
    while( abs((COMs[-1] - COMs[-2])/COMs[-1]) > 10**(-4)*s_diameter):
        for k in range(100):
                        for p in sphere_overlap(distance_balls(balls), balls):
                            if sm_s[str(p)] > 10**(-2)*s_diameter:
                                Cg -= 1.25*10**(-7)
                                if sphere_bottom_overlap(balls)[str(p)] < 0:
                                    balls[str(p)].position_change_gravity(Cg)
                                    balls[str(p)].cart_to_cylindrical()




                                if sphere_wall_overlap(balls)[str(p)] > 0:
                                    balls[str(p)].position_change_cylinder( sphere_wall_overlap(balls), p )
                                    balls[str(p)].cart_to_cylindrical()

                                if sphere_bottom_overlap(balls)[str(p)] > 0:
                                        balls[str(p)].position_change_bottom( sphere_bottom_overlap(balls), p )
                                        balls[str(p)].cart_to_cylindrical()


                                for u in sphere_overlap(distance_balls(balls),balls)[str(p)]:
                                        if sphere_overlap(distance_balls(balls),balls)[str(p)][str(u)] > 0:
                                            balls[str(p)].position_change(unit_vectors(balls, distance_balls(balls)), sphere_overlap(distance_balls(balls),balls),p,u)
                                            balls[str(p)].cart_to_cylindrical()


                                sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
                                sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
                                sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]




                            if sm_s[str(p)] < 0.5*10**(-2)*s_diameter:
                                Cg += 1.25*10**(-7)
                                if sphere_bottom_overlap(balls)[str(p)] < 0: 
                                    balls[str(p)].position_change_gravity(Cg)
                                    balls[str(p)].cart_to_cylindrical()


                                if sphere_wall_overlap(balls)[str(p)] > 0:
                                    balls[str(p)].position_change_cylinder( sphere_wall_overlap(balls), p )
                                    balls[str(p)].cart_to_cylindrical()

                                if sphere_bottom_overlap(balls)[str(p)] > 0:
                                        balls[str(p)].position_change_bottom( sphere_bottom_overlap(balls), p )
                                        balls[str(p)].cart_to_cylindrical()


                                for u in sphere_overlap(distance_balls(balls),balls)[str(p)]:
                                        if sphere_overlap(distance_balls(balls),balls)[str(p)][str(u)] > 0:
                                            balls[str(p)].position_change(unit_vectors(balls, distance_balls(balls)), sphere_overlap(distance_balls(balls),balls),p,u)
                                            balls[str(p)].cart_to_cylindrical()


                                sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
                                sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
                                sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]



                            #print(sm_c)
                            #print(sm_s)
                            #print(sm_b)
                 #       int_dict[p+1] = centre_of_mass(balls)
        COMs = np.append(COMs, centre_of_mass_length(centre_of_mass(balls)))
        print(Cg, COMs)
    
    return(Cg, COMs)
    
    
    
    
def final_phase(balls, Cg):
    
    s_diameter = 2*balls["1"].s_radius


    sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
    sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
    sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]
    
    
    COMs = np.array([])
    
    
    
    C = Cg
    
    while (max(list(sm_s.values())) > s_diameter*10**(-4)) or max(sm_c) > s_diameter*10**(-4) or max(sm_b) > s_diameter*10**(-4):
        COMs = np.append(COMs, 0)
        COMs = np.append(COMs, centre_of_mass_length(centre_of_mass(balls)))
    
        while( abs((COMs[-1] - COMs[-2])/COMs[-1]) > s_diameter*10**(-4)):
        
            for k in range(100):
                            for p in sphere_overlap(distance_balls(balls), balls):
                                    if sphere_bottom_overlap(balls)[str(p)] < 0:
                                        balls[str(p)].position_change_gravity(C)
                                        balls[str(p)].cart_to_cylindrical()




                                    if sphere_wall_overlap(balls)[str(p)] > 0:
                                        balls[str(p)].position_change_cylinder( sphere_wall_overlap(balls), p )
                                        balls[str(p)].cart_to_cylindrical()

                                    if sphere_bottom_overlap(balls)[str(p)] > 0:
                                            balls[str(p)].position_change_bottom( sphere_bottom_overlap(balls), p )
                                            balls[str(p)].cart_to_cylindrical()


                                    for u in sphere_overlap(distance_balls(balls),balls)[str(p)]:
                                            if sphere_overlap(distance_balls(balls),balls)[str(p)][str(u)] > 0:
                                                balls[str(p)].position_change(unit_vectors(balls, distance_balls(balls)), sphere_overlap(distance_balls(balls),balls),p,u)
                                                balls[str(p)].cart_to_cylindrical()


                                    sm_c = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["cylinder"]
                                    sm_s = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["sphere"]
                                    sm_b = sigma_max(sphere_overlap(distance_balls(balls), balls), sphere_wall_overlap(balls), sphere_bottom_overlap(balls))["bottom"]


            COMs = np.append(COMs, centre_of_mass_length(centre_of_mass(balls)))
            print(C, COMs, abs((COMs[-1] - COMs[-2])/COMs[-1]), sm_c, max(list(sm_s.values())), sm_b)
        C = C/2.0
    
    return(C, COMs)
       
    





def bridge_s_s(balls):
    b_s_s = {}
    for p in sphere_overlap(distance_balls(balls),balls):
        for u in sphere_overlap(distance_balls(balls),balls)[str(p)]:
            if sphere_overlap(distance_balls(balls),balls)[str(p)][str(u)] > -2*10**(-2):
                b_s_s[p + "-" + u] = {"x" : { str(p) : balls[str(p)].position_cartesian["x"] , str(u) : balls[str(u)].position_cartesian["x"] },
                                      "y" : { str(p) : balls[str(p)].position_cartesian["y"] , str(u) : balls[str(u)].position_cartesian["y"] },
                                      "z" : { str(p) : balls[str(p)].position_cartesian["z"] , str(u) : balls[str(u)].position_cartesian["z"] }}

                
    ll = list(b_s_s.keys()) 
    n_l = []
    
    for u in list(b_s_s.keys()):
        ll.remove(u)

        for p in ll:
            if (sum(list(b_s_s[u]['x'].values())) ==  sum(list(b_s_s[p]['x'].values()))):
                n_l.append(u)
                
    for k in n_l:
        del b_s_s[k]

        

    
    return(b_s_s)






def bridge_s_w(balls):
    s_radius = balls["1"].s_radius
    b_s_w = {}
    for p in sphere_overlap(distance_balls(balls),balls):
        if sphere_wall_overlap(balls)[str(p)] > 0.0:
            b_s_w[p] = { str(p) : {"x" : balls[str(p)].position_cartesian["x"],
                        "y" : balls[str(p)].position_cartesian["y"] ,
                        "z" : balls[str(p)].position_cartesian["z"] },
                        "walls" + str(p) : { "x" : balls[str(p)].position_cartesian["x"] + 1.1*s_radius*np.cos(balls[str(p)].position["phi"]),
                                             "y" : balls[str(p)].position_cartesian["y"] + 1.1*s_radius*np.sin(balls[str(p)].position["phi"]),
                                             "z" : balls[str(p)].position_cartesian["z"]}}
                        
    return(b_s_w)
        
        
        
def bridge_s_b(balls):
    s_radius = balls["1"].s_radius
    b_s_b = {}
    for p in sphere_overlap(distance_balls(balls),balls):
        if sphere_bottom_overlap(balls)[str(p)] > -2*10**(-4):
            b_s_b[p] = { str(p) : {"x" : balls[str(p)].position_cartesian["x"],
                        "y" : balls[str(p)].position_cartesian["y"] ,
                        "z" : balls[str(p)].position_cartesian["z"] },
                        "bottom" + str(p) : { "x" : balls[str(p)].position_cartesian["x"] ,
                                             "y" : balls[str(p)].position_cartesian["y"],
                                             "z" : balls[str(p)].position_cartesian["z"] - 1.1*s_radius}}
                        
    return(b_s_b)
        
