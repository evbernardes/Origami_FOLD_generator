#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:29:48 2020

@author: evandro
"""
import numpy as np
from .Origami import Origami
# from Origami import Origami

# pi,sin,cos,tan,arcsin,arccos,arctan,arctan2 = np.pi,np.sin,np.cos,np.tan,np.arcsin,np.arccos,np.arctan,np.arctan2
TORAD = np.pi/180.0
TODEG = 1/TORAD

class Bendy(Origami):
    def __init__(self,R,ratio,n,m,alpha1,alpha2,top_closed = False,bottom_closed = False,cylinder_closed=True,bent = False,base_height = 0):
        
        Origami.__init__(self,frame_title = 'Bendy')
        self.frame_attributes.append('nonSelfTouching')
        
        alpha1 = alpha1*TORAD
        alpha2 = alpha2*TORAD
        beta1 = np.arccos(np.cos(alpha1)*np.sin(np.pi/n))
        beta2 = np.arccos(np.cos(alpha2)*np.sin(np.pi/n))
        b1 = R*(1-ratio)*np.sin(beta1)/np.cos(alpha1)
        b2 = R*(1-ratio)*np.sin(beta2)/np.cos(alpha2)
        apo = R*np.cos(np.pi/n)
        A = 2*R*np.sin(np.pi/n)
        a = ratio*A
        R_ = np.sqrt(a*a/4 + apo*apo)
        phi = np.pi/n - abs(np.arctan((a/2)/apo))

        ring_stages = []
        middle_stages_1 = []
        middle_stages_2 = []
        
        # m = 1
        for j in range(m+1):
            ring_stage = []
            for i in range(n):
                coords = [R*np.cos(i*2*np.pi/n),R*np.sin(i*2*np.pi/n),(b1+b2)*j + base_height]
                ring_stage.append('ring_'+str(j)+'_'+str(i))
                self.add_vertex(ring_stage[-1],coords)
            ring_stages.append(ring_stage)
            
            if j < m:  
                middle_stage_1 = []
                middle_stage_2 = []
                for i in range(n):
                    r = R_ if bent and i == 0 else R
                    # r = R
                    coords = [r*np.cos(i*2*np.pi/n),r*np.sin(i*2*np.pi/n),(b1+b2)*j + b1 + base_height]
                    middle_stage_1.append('middle1_'+str(j)+'_'+str(i))
                    self.add_vertex(middle_stage_1[-1],coords)
                    r = R if (bent and i == 0) else R_
                    # r = R_
                    coords = [r*np.cos(i*2*np.pi/n + phi),r*np.sin(i*2*np.pi/n + phi),(b1+b2)*j + b1 + base_height]
                    middle_stage_2.append('middle2_'+str(j)+'_'+str(i))
                    self.add_vertex(middle_stage_2[-1],coords)
                middle_stages_1.append(middle_stage_1)
                middle_stages_2.append(middle_stage_2)
        
        # top and bottom ring mountain folds
        for i in range(n):
            bottom_type = 'M' if bottom_closed or base_height>0 else 'B'
            self.connect(ring_stages[0][i],ring_stages[0][(i+1)%n],bottom_type)
            top_type = 'M' if top_closed  or base_height>0 else 'B'
            self.connect(ring_stages[-1][i],ring_stages[-1][(i+1)%n],top_type)
        
        for j in range(m):
            for i in range(n):
                # middle mountain fold rings
                if j < m-1:
                    self.connect(ring_stages[j+1][i],ring_stages[j+1][(i+1)%n],'M')
                    
                self.connect(middle_stages_1[j][i],middle_stages_2[j][i],'V' if bent and i == 0 else 'M')
                # self.connect(middle_stages_2[j][i],middle_stages_1[j][(i+1)%n],'M' if bent and i == 0 else 'V')
                self.connect(middle_stages_2[j][i],middle_stages_1[j][(i+1)%n],'F' if bent and i == 0 else 'V')
                self.connect(ring_stages[j][i],middle_stages_2[j][(i)%n],'M' if bent and i == 0 else 'V')
                self.connect(ring_stages[j+1][i],middle_stages_2[j][(i)%n],'M' if bent and i == 0 else 'V')
                
                if cylinder_closed == False and i == n-1:
                    self.connect(ring_stages[j][i],middle_stages_1[j][i],'B')
                    self.connect(middle_stages_1[j][i],ring_stages[j+1][i],'B')
                elif bent == True and i == 0:
                    # self.connect(ring_stages[j][i],middle_stages_1[j][i],'F')
                    # self.connect(middle_stages_1[j][i],ring_stages[j+1][i],'F')
                    self.connect(ring_stages[j][i],middle_stages_1[j][i],'U')
                    self.connect(middle_stages_1[j][i],ring_stages[j+1][i],'U')
                else:
                    self.connect(ring_stages[j][i],middle_stages_1[j][i],'M')
                    self.connect(middle_stages_1[j][i],ring_stages[j+1][i],'M')
                
                if True:
                    face_left_bottom = [ring_stages[j][i],middle_stages_2[j][i],middle_stages_1[j][i]]
                    self.add_face(face_left_bottom)
                    face_left_top = [middle_stages_1[j][i],middle_stages_2[j][i],ring_stages[j+1][i]]
                    self.add_face(face_left_top)
                    face_right_bottom = [ring_stages[j][i],ring_stages[j][(i+1)%n],middle_stages_1[j][(i+1)%n],middle_stages_2[j][(i+0)%n]]
                    self.add_face(face_right_bottom)
                    face_right_top = [middle_stages_2[j][(i+0)%n],middle_stages_1[j][(i+1)%n],ring_stages[j+1][(i+1)%n],ring_stages[j+1][i]]
                    self.add_face(face_right_top)
          
        if base_height > 0:
            base_top = []
            base_bottom = []
            for i in range(n):
                coords = [R*np.cos(i*2*np.pi/n),R*np.sin(i*2*np.pi/n),0]
                base_bottom.append('base_bottom_'+str(i))
                self.add_vertex(base_bottom[i],coords)
                self.connect(base_bottom[i],ring_stages[0][i],'F',angle = np.pi*(n-2)/n)
                
                coords[2] = coords[2] + base_height*2 + m*(b1+b2)
                base_top.append('base_top_'+str(i))
                self.add_vertex(base_top[i],coords)
                self.connect(base_top[i],ring_stages[-1][i],'F',angle = np.pi*(n-2)/n)
                
            top_color = 'M' if top_closed else 'B'
            bottom_color = 'M' if bottom_closed else 'B'
            for i in range(n):
                self.connect(base_top[i],base_top[(i+1)%n],top_color)
                self.connect(base_bottom[i],base_bottom[(i+1)%n],bottom_color)
                
                face_base_top = [ring_stages[-1][i],ring_stages[-1][(i+1)%n],base_top[(i+1)%n],base_top[i]]
                self.add_face(face_base_top)
                face_base_bottom = [base_bottom[i],base_bottom[(i+1)%n],ring_stages[0][(i+1)%n],ring_stages[0][i]]
                self.add_face(face_base_bottom)
        else:
            base_top = ring_stages[-1]
            base_bottom = ring_stages[0]
                    
        # top and bottom facets
        if top_closed:
            self.add_vertex('top_center',[0,0,base_height*2 + m*(b1+b2)])
            for i in range(n):
                self.connect('top_center',base_top[i],'F')
                self.add_face(['top_center',base_top[i],base_top[(i+1)%n]])
        if bottom_closed:
            self.add_vertex('bottom_center',[0.,0.,0.])
            for i in range(n):
                self.connect('bottom_center',base_bottom[i],'F')
                self.add_face(['bottom_center',base_bottom[(i+1)%n],base_bottom[i]]) 


# bendy = Bendy(1,0.4,8,2,35,45).plot(faces=True)     
# bendy = Bendy(40,0.4,8,4,35,45,bent = True,top_closed=True,bottom_closed=True,base_height = 0)
# bendy.plot(faces=False)           
# bendy.save("/home/evandro/Desktop/fold_test/bendy_numpy.fold")