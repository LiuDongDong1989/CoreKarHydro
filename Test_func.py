# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 14:32:01 2021

@author: liudongdong

"""
import numpy as np

def test_func(parameters, treatment=None):
    a, b,=parameters
    x= np.array([1,2,3,4])
    y=a*x+b
    return y







