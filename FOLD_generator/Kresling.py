#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:29:48 2020

@author: evandro
"""
import numpy as np
from .Origami import Origami

class Kresling(Origami):
    def __init__(self,R,n,m,lamb = 0.7,top_closed = False,bottom_closed = False,cylinder_closed=True):
        
        Origami.__init__(self,frame_title = 'Kresling')
        self.frame_attributes.append('nonSelfTouching')
        
        theta = np.pi*(n-2)/(2*n)
        a = 2*R*np.sin(np.pi/n)
        l = 2*R*np.cos(theta*(1 - lamb))
        b = np.sqrt(l**2 + a**2 - 2*l*a*np.cos(theta*lamb))
        h = np.sqrt((l*l + b*b) - 4 * R*R)
        alpha = (np.arccos(1 - 2*(l*l - h*h)/(4*R*R)) - 4*np.pi/n)
        
        stages = []
        
        for j in range(m+1):
            stage = []
            for i in range(n):
                coords = (R*np.cos(i*2*np.pi/n+alpha*j),R*np.sin(i*2*np.pi/n+alpha*j),h*j)
                stage.append('v_'+str(j)+'_'+str(i))
                self.add_vertex(stage[-1],coords)
            stages.append(stage)
        
        # top and bottom ring mountain folds
        for i in range(n):
            bottom_type = 'M' if bottom_closed else 'B'
            self.connect(stages[0][i],stages[0][(i+1)%n],bottom_type)
            top_type = 'M' if top_closed else 'B'
            self.connect(stages[-1][i],stages[-1][(i+1)%n],top_type)
        
        # top and bottom facets
        if top_closed:
            self.add_face(stages[-1])
        if bottom_closed:
            self.add_face(stages[0][::-1])
        
        for j in range(m):
            for i in range(n):
                # middle mountain fold rings
                if j < m-1:
                    self.connect(stages[j+1][i],stages[j+1][(i+1)%n],'M')
          
                # mountain and valley diagonal folds
                self.connect(stages[j][(i-2)%n],stages[j+1][i],'V')
                if cylinder_closed == False and i == n-1:
                    self.connect(stages[j][(i-1)%n],stages[j+1][i],'B')
                else:
                    self.connect(stages[j][(i-1)%n],stages[j+1][i],'M')
                
                # every triangular facet
                facet_a = [stages[j+1][(i+1)%n],stages[j+1][i],stages[j][(i-1)%n]]
                facet_b = [stages[j+1][(i+1)%n],stages[j][(i-1)%n],stages[j][(i)%n]]
                self.add_face(facet_a)
                self.add_face(facet_b)
                
        
# def __main__():
#     Kresling(1,8,1,0.7).plot()