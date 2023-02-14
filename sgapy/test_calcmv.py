# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:48:05 2023

@author: bmelosh
"""
import numpy as np

from calcmv import calcmv

T = [1,2,3]
P = [1,2,3]

print(np.degrees(T))
print(np.degrees(P))

calcmv(T,P)
