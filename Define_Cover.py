# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:27:47 2019
#---------------------函数说明-------------------------------------------------#
#定义覆盖区域
#-----------------------------------------------------------------------------#
@author: Lenovo
"""
def Define_Cover(mesh, X, Y, Lx, Ly, coverage=0): 
    BSCArea=Lx*Ly*coverage#苔藓区域的面积[m]
    BSCCover=[Lx/2-pow(BSCArea,0.5)/2, Lx/2+pow(BSCArea,0.5)/2]#苔藓的坐标
    mask=( (X< BSCCover[0]) | (X > BSCCover[1]) )#定义特殊区域
    Soilfaces=mesh.facesTop & mask#土壤区域
    BSCfaces= mesh.facesTop & ~mask#苔藓覆盖区域 
    return Soilfaces, BSCfaces