# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:27:50 2019

@author: Lenovo
"""
#---------------------函数说明-------------------------------------------------#
#自定义大气边界条件
#上边界蒸发速率随表土含水量的改变发生变化
#psi:土壤水吸力
#psi_hygro:土壤吸湿点
#psi_field:土壤田持含水量下的水吸力
#-----------------------------------------------------------------------------#
import numpy as np 

def UpBoundary(mesh, X, Y, timeStepDuration,
              psi, theta, Krel, psi_tmp, psi_max, 
              PET, P, AET_Soil, AET_BSC,
              PNfactor, BSCfactor,
              psi_hygro, psi_field,            
              Soilfaces, BSCfaces,
              step=0):
    #定义蒸发速率函数公式的参数
    AET_max=PET[step]
    
    #土壤区的大气边界条件(仅被松针覆盖)
    for i in range(np.where(Soilfaces.value==1)[0].size):
        point=np.where(Soilfaces.value==1)[0][i]
        criticalValue=psi.faceValue.value[point]
        mask= ((X==X[point]) & (Y==Y[point]))
        subRegion1=(Soilfaces & mask)
        # print(i, "土壤区的", criticalValue)
        if criticalValue<=psi_hygro:
           #在外界的蒸发力超过土壤的输水能力时, 土表处于风干状态,此时的边界条件为水头边界
           psi.constrain(psi_hygro, subRegion1) 
        elif (psi_hygro<criticalValue<=0):
           #当灌溉/降雨强度 ,或蒸发强度并未超出土壤入渗/蒸发能力时, 上边界条件为流量边界条件 
           if (psi_hygro<criticalValue<=psi_field):  
               #表层含水量处于吸湿点和田持点之间时，土壤蒸发速率与含水率成线性关系
               AET_Soil.faceValue[point]=AET_max*PNfactor*(criticalValue-psi_hygro)/(psi_field-psi_hygro)          
           elif (psi_field<criticalValue<=0):
               #表层含水量大于田持点之间时，土壤蒸发速率达到最大 
               AET_Soil.faceValue[point]=AET_max*PNfactor
           flux=(AET_Soil.faceValue[point]-P[step]+Krel.faceValue[point]*[[0],[1]])/Krel.faceValue[point]*(-1)
           psi.faceGrad.constrain(flux, subRegion1)         
        elif (0<criticalValue<=psi_max):
           #当降雨\灌溉强度超出土壤入渗能力时 ,地表将会出现积水现象, 积水深度在没有超出最大深度   
           diffpsi=(psi.faceValue[point]-psi_tmp.faceValue[point])/timeStepDuration
           AET_Soil.faceValue[point]=0
           flux=(AET_Soil.faceValue[point]-P[step]+Krel.faceValue[point]*[[0],[1]]+diffpsi)/Krel.faceValue[point]*(-1)
           psi.faceGrad.constrain(flux, subRegion1)   
        else:
           #当积水深度超出地表最大滞蓄深度时 ,积水深度不再增加, 超渗的水量即形成地表径流
           psi.constrain(psi_max, subRegion1)   
           
    #覆盖区的大气边界条件（被苔藓和松针同时覆盖）
    for i in range(np.where(BSCfaces.value==1)[0].size):
        point=np.where(BSCfaces.value==1)[0][i]
        criticalValue=psi.faceValue.value[point]
        mask= ((X==X[point]) & (Y==Y[point]))
        subRegion2=(BSCfaces & mask)
        if criticalValue<=psi_hygro:
           #在外界的蒸发力超过土壤的输水能力时, 土表处于风干状态,此时的边界条件为水头边界
           psi.constrain(psi_hygro, subRegion2) 
        elif (psi_hygro<criticalValue<=0): 
           #当灌溉/降雨强度 ,或蒸发强度并未超出土壤入渗/蒸发能力时, 上边界条件为流量边界条件 
           if (psi_hygro<criticalValue<=psi_field):  
               #表层含水量处于吸湿点和田持点之间时，土壤蒸发速率与含水率成线性关系  
               AET_BSC.faceValue[point]=AET_max*PNfactor*BSCfactor*(criticalValue-psi_hygro)/(psi_field-psi_hygro)
           elif (psi_field<criticalValue<=0):
               #表层含水量大于田持点之间时，土壤蒸发速率达到最大     
               AET_BSC.faceValue[point]=AET_max*PNfactor*BSCfactor
           flux=(AET_BSC.faceValue[point]-P[step]+Krel.faceValue[point]*[[0],[1]])/Krel.faceValue[point]*(-1)
           psi.faceGrad.constrain(flux, subRegion2)          
        elif (0<criticalValue<=psi_max):
           #当降雨\灌溉强度超出土壤入渗能力时 ,地表将会出现积水现象, 积水深度在没有超出最大深度   
           diffpsi=(psi.faceValue[point]-psi_tmp.faceValue[point])/timeStepDuration 
           AET_BSC.faceValue[point]=0
           flux=(AET_BSC.faceValue[point]-P[step]+Krel.faceValue[point]*[[0],[1]]+diffpsi)/Krel.faceValue[point]*(-1)
           psi.faceGrad.constrain(flux, subRegion2)           
        else:
           #当积水深度超出地表最大滞蓄深度时 ,积水深度不再增加, 超渗的水量即形成地表径流
           psi.constrain(psi_max, subRegion2) 
    return psi  