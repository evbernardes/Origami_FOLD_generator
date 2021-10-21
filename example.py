#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 12:21:45 2020

@author: evandro
"""
import FOLD_generator as fold

R = 25.0/2
ratio = 0.4
n = 8
m = 1
alpha1 = 45.
alpha2 = alpha1 - 10.
base_height = 0

origami = fold.Bendy(R,ratio,n,m,alpha1,alpha2,bent = False,top_closed=False,bottom_closed=False,base_height = base_height)
origami.plot(faces=False)
origami.save('test.fold')


