# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 17:11:37 2019
#---------------------函数说明-------------------------------------------------#
#输入气象数据等实测数据
#-----------------------------------------------------------------------------#
@author: Lenovo
"""
from pandas import read_excel
import numpy as np

def Meteo(Lx, Ly, treatment='B1P1'):
    #时间序列的气象数据等实测数据
    # meteo=read_excel('InputData.xls', header=0, sep = '\t', encoding = 'utf-8')
    meteo=read_excel('InputData.xls', header=0)
    Time=meteo["Time(s)"]#时间
    Tair=meteo["Tair"]#大气温度
    RH=meteo["RH"]#大气湿度
    TopSurface=np.pi*(Lx/2)**2
    UptakezoneVolume=np.pi*(Lx/2)**2*Ly
    P=meteo["P"]/1e6/(3600*24)/UptakezoneVolume#降雨强度[1/s]
    PET=meteo["PET"]/1000/(3600*24)*TopSurface/UptakezoneVolume#实测大气蒸发速率[1/s]
    SWC_true=meteo[treatment]#实测土壤平均含水率[-]
    return [Time, Tair, RH, P, PET, SWC_true]

"""
Created on Wed Dec  4 08:03:29 2019
#---------------------函数说明-------------------------------------------------#
#计算潜在蒸发速率PET、覆盖影响下的蒸发计算
# Actual Evapotranspiration [s-1]. Computed by the solver.
# Potential Evapotranspiration [s-1].in each cell of the uptake zone :
#(1)PET = measured value of PET [m.s-1] * Top surface [m2] / Uptake zone volume [m3]
#(2)anywhere else (i.e., deeper than uptake zone) : PET = 0.
#-----------------------------------------------------------------------------#
@author: Lenovo
"""
def Calculation_PET(Tair=0, RH=0, Lx=0):        
    es=0.6108*10**(7.5*Tair/(273.3+Tair))
    ea=es*RH/100
    PET_pre=0.85*(es-ea)/1000/(3600*24*Lx)
    return PET_pre