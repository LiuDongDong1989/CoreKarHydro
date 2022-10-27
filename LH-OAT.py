# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 14:45:26 2021
目的：LH-OAT敏感性分析方法
编程思路：
模块1：拉丁超立方抽样
模块2：某个LH抽样点执行OAT算法
模块3：全局敏感性分析
@author: Lenovo
"""
import MOCOM_UA_class
from SALib.sample import latin

class LH_OAT(MOCOM_UA_class.ParetoSet): 
    
     def __init__(self,problem):    
         MOCOM_UA_class.ParetoSet.__init__(self, problem)#继承父类
         self.ParameterSpace= latin.sample(problem, 500)#LH拉丁超立方抽样  
         
     #某个LH抽样点执行OAT算法    
     def OAT(self):
         self.fraction=0.05
         for i, element in enumerate(self.ParameterSpace):
             self.ObjectiveSpace[i, :], self.ParameterSpace[i] =self.cal_objectives_single(element)  
             M_change
         S_ij=abs( ( (M_change-M)/((M_change-M)/2) ) / self.fraction)
         
         
         
         
    
    
    
    