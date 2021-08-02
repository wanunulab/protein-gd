# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 14:40:28 2020

@author: Ali
"""
import numpy as np
def load_bdat(filename,msec_steps=0.004):
    current=np.fromfile(filename) #* 1000
    return msec_steps,current 