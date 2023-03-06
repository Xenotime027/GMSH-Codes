#This script creates a 2D material file with GMSH format containing a homogenous matrix with a from user set number of inclusions
#This is useful for automatic generation of finite-element files and can be modified to include inclusions with other shapes and morphology.

import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import datetime

def generate_gmsh(iter, num_min=2, num_max=9, radius_mu=0.3, radius_sigma=0.03, box_size=[3,3], num_mesh_circle=51, num_mesh_box=31):
    #inputs:
        #iter: the current number of iteration
        #num_min: minimum number of inclusions 
        #num_max: maximum number of inclusions
        #radius_mu: mean radius for the circles
        #radius_sigma: standard deviation for the radius
        #box_size: dimensions of the box (can be square or rectangular, preferable to use positive values)
        #num_mesh_circle: number of points on each circle to place nodal points
        #num_mesh_box: number of nodal points on each side of the box
        
    #num: generates a random number of inclusion within the range of num_min and num_max
    num = random.randint(num_min, num_max)
    #counter: counter until the number of inclusions is reached
    counter = num
    #trial: counter to force end the circle generation if overlap threshold is reached (to save time) 
    trial = 0
    #num_circle: counter of the current number of circles
    num_circle = 0
    #collector: to collect the coordinates (x and y) of the centre of each circle
    collector = np.zeros((num, 2))
    #radius: a normal distribution of radii, calculated from radius_mu, radius_sigma
    radius = np.random.normal(radius_mu, radius_sigma, num)

    while counter:
        ob = True
        coord = []
        coord.append(np.random.uniform(1.02*radius[num_circle], box_size[0]-1.02*radius[num_circle]))
        coord.append(np.random.uniform(1.02*radius[num_circle], box_size[1]-1.02*radius[num_circle]))
        
        #force end mechanism if no more circles can be found (due to overlapping)
        trial +=1
        if trial > 3333:
            ob = False
            break
        
        #forces the circle to be inside the box
        for i in range(2):
            if(coord[i]-2*radius[num_circle]<0 or coord[i]+2*radius[num_circle]>box_size[i]):
                ob = False
                break
        
        #forces the circles not to overlap with one another with euclidean distance
        for j in range(num_circle):
            dist_x = distance.euclidean(coord, collector[j])  
            
            if dist_x < (3*(radius[j]+radius_sigma)):
                ob = False
                break
        #if everything is successful, the centre coordinates are added to the collector
        if ob:
            counter -= 1   
            collector[num_circle,:] = coord
            num_circle += 1
    
    #uncomment this to show the circles in plot (do not use if the number of iterations is large)
    fig, ax = plt.subplots(figsize=(2*box_size[0], 2*box_size[1]))
    ax.set_xlim((0, box_size[0]))
    ax.set_ylim((0, box_size[1]))

    for i in range(num_circle):
        circle = plt.Circle((collector[i][0], collector[i][1]), radius[i])
        ax.add_patch(circle)

    plt.show()
            
    #writes the input file
    with open(f"{iter}.geo", "w") as f: 
        #print(iter)
        dt = datetime.datetime.today()
        f.write(f'// Gmsh project created on {dt.strftime("%a")} {dt.strftime("%b")} {dt.day} {dt.hour}:{dt.minute}:{dt.second} {dt.year}')
        f.write('\nSetFactory("OpenCASCADE");')
        f.write(f"\n//+")
        
        #defines the circles
        for i in range(num_circle):
            f.write(f"\nCircle({i+1}) = {{{collector[i][0]}, {collector[i][1]}, 0, {radius[i]}, 0, 2*Pi}};")
            f.write(f"\n//+")
        
        point = num_circle+1
        line = num_circle+1
        
        #defines the points of the box
        f.write(f"\nPoint({point}) = {{0, 0, 0, 1.0}};\n//+")
        f.write(f"\nPoint({point+1}) = {{0, {box_size[1]}, 0, 1.0}};\n//+")
        f.write(f"\nPoint({point+2}) = {{{box_size[0]}, 0, 0, 1.0}};\n//+")
        f.write(f"\nPoint({point+3}) = {{{box_size[0]}, {box_size[1]}, 0, 1.0}};\n//+")
        
        #defines the lines of the box
        f.write(f"\nLine({line}) = {{{point+2}, {point}}};\n//+")
        f.write(f"\nLine({line+1}) = {{{point}, {point+1}}};\n//+")
        f.write(f"\nLine({line+2}) = {{{point+1}, {point+3}}};\n//+")
        f.write(f"\nLine({line+3}) = {{{point+3}, {point+2}}};\n//+")
        
        #defines the surfaces of the inclusions
        for i in range(num_circle):
            f.write(f"\nCurve Loop({i+1}) = {{{i+1}}};")
            f.write(f"\n//+")
            f.write(f"\nPlane Surface({i+1}) = {{{i+1}}};")
            f.write(f"\n//+")
        
        curve_loop = num_circle + 1
        
        #defines the surface of the box
        f.write(f"\nCurve Loop({curve_loop}) = {{{line+2}, {line+3}, {line}, {line+1}}};")
        f.write(f"\n//+")
        
        #defines the excluded surface of the box (not to include the inclusions)
        for i in range(num_circle):
            f.write(f"\nCurve Loop({i+curve_loop+1}) = {{{i+1}}};")
            f.write(f"\n//+")
        
        #string generator to exclude the inclusions
        x = str(curve_loop)+str(", ")
        for i in range(num_circle):
            x = x + str(curve_loop+i+1)+str(", ")
        x = x[:-2]
        
        f.write(f"\nPlane Surface({curve_loop}) = {{{x}}};")
        f.write(f"\n//+")
        
        #defines the physical curves and their names
        p_curve_name = ["bottom", "left", "top", "right"]
        p_curve = curve_loop + num_circle + 1
        
        #box edges as physical curves
        for i in range(4):
            f.write(f'\nPhysical Curve("{p_curve_name[i]}", {p_curve+i}) = {{{line+i}}};')
            f.write(f"\n//+")
        
        for i in range(num_circle):
            f.write(f"\nTransfinite Curve {{{i+1}}} = {num_mesh_circle} Using Progression 1;")
            f.write(f"\n//+")

        for i in range(num_circle+1, num_circle+5):
            f.write(f"\nTransfinite Curve {{{i}}} = {num_mesh_box} Using Progression 1;")
            f.write(f"\n//+")
        
        #string generator to define physical surface for the inclusion
        y = str()
        for i in range(num_circle):
            y = y + str(i+1) + str(", ")
        y = y[:-2]
        
        #defines the physical surfaces of the inclusion and matrix
        f.write(f'\nPhysical Surface("inclusion", {p_curve+4}) = {{{y}}};')
        f.write(f"\n//+")
        f.write(f'\nPhysical Surface("matrix", {p_curve+5}) = {curve_loop};')
        f.write(f"\n")
        
        
        
        
###############################################
#Example implementation
#Example Parameters for generating inclusions
num_min = 11
num_max = 12
radius_mu = 0.25
radius_sigma = radius_mu/10
box_size = [5, 5]
num_mesh_circle = 21
num_mesh_box = 31
num_iter = 1

#Implementation
for iter in range(num_iter):
    generate_gmsh(iter, 
                  num_min=num_min, 
                  num_max=num_max, 
                  radius_mu=radius_mu,
                  radius_sigma=radius_sigma, 
                  box_size=box_size, 
                  num_mesh_circle=num_mesh_circle, 
                  num_mesh_box=num_mesh_box)
    
