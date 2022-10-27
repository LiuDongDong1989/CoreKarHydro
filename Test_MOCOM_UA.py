# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:43:35 2021
用于测试参数反演程序
@author: liudongdong
"""
#测试MOCOM_UA代码 
import os
os.environ.setdefault("FIPY_SOLVERS", "scipy")   
import MOCOM_UA_class

'''
#problem说明：
#:param num_vars:参数个数
#:param names:参数名称（自定义）
#:param bounds:参数对应范围（list）
#:param objNum:多目标函数个数
#:param population:种群数量
#:param func:正演程序 
''' 

flag='VadoseZone'
if flag=='Test_func':
    import Test_func #常规线性拟合测试
    problem = {'num_vars': 2,
               'names': ['x0', 'x1'],
               'bounds': [[1.5, 2.5], [0.5, 1.5]],
               'objNum': 2,
               'population': 100,
               'func': Test_func.test_func,
               'outputs': flag,
               'treatment': None
               }
elif flag=='VadoseZone':
    import VadoseZone #VadoseZone测试等通用模型
    treatment='B0P1'
    if treatment[1]=='0' and treatment[-1]=='0':
        problem= {'num_vars': 5,
           'names': ['thetas', 'thetar','alpha', 'n', 'Ks'],
           'bounds': [[0.4, 0.6], [0.0, 0.2],[1, 6], [1.1, 1.7],[1e-8, 1e-4]],
           'objNum': 2,
           'population': 100,
           'func': VadoseZone.VadoseZone,
           'outputs': flag,
           'treatment': treatment} 
    else:
        if treatment[1]!='0' and treatment[-1]!='0':
            problem= {'num_vars': 4,
                       'names': ['BSCfactor', 'PNfactor','DecayFactor', 'DecayRate'],
                       'bounds': [[0.3, 1], [0.3, 1], [0.3, 1], [0, 0.5]],
                       'objNum': 2,
                       'population': 100,
                       'func': VadoseZone.VadoseZone,
                       'outputs': flag,
                       'treatment': treatment}
        elif treatment[1]!='0' and treatment[-1]=='0':
            problem= {'num_vars': 3,
                       'names': ['BSCfactor','DecayFactor', 'DecayRate'],
                       'bounds': [[0.3, 1], [0.3, 1], [0, 0.5]],
                       'objNum': 2,
                       'population': 100,
                       'func': VadoseZone.VadoseZone,
                       'outputs': flag,
                       'treatment': treatment}
        elif treatment[1]=='0' and treatment[-1]!='0':
            problem= {'num_vars': 1,
                       'names': ['PNfactor'],
                       'bounds': [[0.3, 1]],
                       'objNum': 2,
                       'population': 100,
                       'func': VadoseZone.VadoseZone,
                       'outputs': flag,
                       'treatment': treatment}
            
optimal=MOCOM_UA_class.ParetoSet(problem)#构建实例对象
optimal.MOCOM_UA() #参数反演

