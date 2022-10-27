# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:42:52 2019

@author: Lenovo
"""
'''---------------------函数说明-------------------------------------------------#
#定义水分特征曲线模型
#    #输入V-G model参数
#    thetas=0.60#saturation water content [-]
#    thetar=0.095#residual water content [-]
#    alpha=1.9#inverse of the capillary length [m-1]
#    n=1.31#van Genuchten exponent [-]
#    S=0.00001#Water storage coefficient [m-1]
#    K=8.844444444444445e-07#Saturated hydraulic conductivity field [m/s]   
#    #方法一：用偏导数的直接公式表示Crel (推荐使用)
#    #方法二：用偏导数的差分公式表示Crel
#    def pos0(a):
#        np.where(a>=0, a, 1)
#        np.where(a<0, a, 0)
#        return a 
#    thtil_tmp=0.5*((1+sign(psi_tmp))+(1-sign(psi_tmp))*pow((1+pow(abs(alpha*psi_tmp),n)),-(1-(1/n))))
#    usf=nx*ny#Dimensionless unit vertical upward vector field[-].
#    Crel=0.5*(1+sign(psi))*S+0.5*((1-sign(psi))*((thetas-thetar)*(thtil-thtil_tmp)*(1./((usf*pos0(psi-psi_tmp)*pos0(psi_tmp-psi))+psi-psi_tmp)))) 
#-----------------------------------------------------------------------------#
'''
from numpy import sign 

def VGmodel(psi, thetas, thetar, alpha, n, Ks, S=1.0E-05):           
    #计算土壤水分特征曲线 V-G模型
    #sign为符号函数
    thtil=0.5*((1+sign(psi))+(1-sign(psi))*pow((1+pow(abs(alpha*psi),n)),-(1-(1/n))))    
    Krel=0.5*((1+sign(psi))*Ks+(1-sign(psi))*Ks*pow(thtil,0.5)*pow((1-pow((1-pow(thtil,(n/(n-1)))),(1-(1/n)))),2))
    Crel=S*0.5*(1+sign(psi))+0.5*((1-sign(psi))*(thetas-thetar)*alpha*n*(1-(1/n))*pow(abs(alpha*psi),(n-1))*pow((1+pow(abs(alpha*psi),n)),-(2-(1/n))))
    theta=(thetas-thetar)*thtil+thetar#上一个时间步长的值 
    # U=(-1)*(Krel.faceValue*psi.faceGrad+Krel.faceValue*[[0],[1]])           
    return [thtil, Krel, Crel, theta]

