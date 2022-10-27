# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:23:34 2019
@author: 刘冬冬 贵州大学资源与环境工程学院
程序说明：正演主程序
"""
import numpy as np
import Define_Cover
import SWRC
import Input_Meteo
import UpBoundary
import math
import pre_Error

def VadoseZone(Parameters, treatment=None, Calibration_percent=0.4):    
    '''
    #参数列表
    #BSCfactor:苔藓覆盖蒸发系数[-]
    #PNfactor:松针覆盖蒸发系数[-]
    #DecayFactor:最终苔藓覆盖系数与初始苔藓覆盖蒸发系数的比值，范围[0,1]
    #DecayRate：退化速率,范围[0,1]
    #psi_field:土壤田持含水量下的水吸力[m]
    #thetas=0.60#saturation water content [-]
    #thetar=0.095#residual water content [-]
    #alpha=1.9#inverse of the capillary length [m-1]
    #n=1.31#van Genuchten exponent [-]
    #S=0.00001#Water storage coefficient [m-1]
    #K=8.844444444444445e-07#Saturated hydraulic conductivity field [m/s]
    #psi_hygro:土壤吸湿点[m]
    #coverage:苔藓覆盖度[-] 
    '''
    print('START...输入参数', Parameters, end='->')
    if treatment[1]=='0' and treatment[-1]=='0':
        BSCfactor_initial=1
        PNfactor=1
        DecayFactor=1
        DecayRate=0
        thetas=Parameters[0]
        thetar=Parameters[1]
        alpha=Parameters[2]
        n=Parameters[3]
        Ks=Parameters[4]             
    else:
        thetas=0.5386
        thetar=0.1393
        alpha=2.611
        n=1.228
        Ks=4.894e-5
        if treatment[1]!='0' and treatment[-1]!='0': #苔藓+松针覆盖
            BSCfactor_initial=Parameters[0]
            PNfactor=Parameters[1]
            DecayFactor=Parameters[2]
            DecayRate=Parameters[3]
        elif treatment[1]!='0' and treatment[-1]=='0': #无松针覆盖
            BSCfactor_initial=Parameters[0]
            PNfactor=1
            DecayFactor=Parameters[1]
            DecayRate=Parameters[2]
        elif treatment[1]=='0' and treatment[-1]!='0': #无苔藓覆盖
            BSCfactor_initial=1
            PNfactor=Parameters[0]
            DecayFactor=1
            DecayRate=0    
    BSCfactor_final=BSCfactor_initial*DecayFactor
    
    import fipy.solvers.solver    
    print('调用求解器(是否采用petsc求解器)...', fipy.solvers.solver,end='->')
   
    print("划分网格...", end='->')
    nx, ny =3, 3
    Lx, Ly = 0.20,0.15
    dx, dy = Lx/nx, Ly/ny #单位：m
    mesh = fipy.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
    timeStepDuration=86400#[s]#timeStepDuration =0.9* dy**2 / (2 * max(Krel))
    X, Y = mesh.faceCenters#单元面的坐标
    
    print("定义覆盖区域、输入气象数据...", end='->')
    coverage={'B0P0':0.00, 'B1P0': 0.1028, 'B2P0':0.2654, 'B3P0': 0.32,
              'B0P1': 0.00,   'B1P1':0.13, 'B2P1': 0.3197, 'B3P1': 0.3402, 
              'B0P2': 0.00,   'B1P2':0.1725, 'B2P2': 0.3385, 'B3P2': 0.3722}
    Time, Tair, RH, P, PET, SWC_true=Input_Meteo.Meteo(Lx, Ly, treatment=treatment)
    Soilfaces, BSCfaces=Define_Cover.Define_Cover(mesh, X, Y, Lx, Ly, coverage=coverage[treatment])
    
    print("初始化水头分布...", end='->')
    # Water pressure field [m]
    psi_max=0.0       #地表最大滞蓄深度[m]
    psi_initial=-10.  #初始水头[m]
    psi=fipy.CellVariable(name="Total water pressure head (m)", mesh=mesh,value=psi_initial, hasOld=1)#当前时间步长的值，扩散系数是非线性的,用hasOld来保持原值直至更新动作
    psi_tmp = fipy.CellVariable(name="Total water pressure head of last time step (m)", mesh=mesh, value=psi_initial)#上一个时间步长的值 
    
    print("定义土壤水力特征参数...", end='->')
    #输入V-G model参数
    psi_hygro=-102/1000*3100 #土壤吸湿点[m]  
    psi_field=-102/1000*30#土壤田持含水量下的水吸力[m]
    thtil, Krel, Crel, theta=SWRC.VGmodel(psi, thetas, thetar,alpha, n, Ks)
      
    print("计算潜水蒸发、初始化蒸发...", end='->')
    # from Calculation_PET import Calculation_PET
    # PET_pre=Calculation_PET(Tair=Tair, RH=RH, Lx=Lx)#大气蒸发量的预测值
    AET_Soil=fipy.CellVariable(name="裸土处的蒸发速率 (s-1)", mesh=mesh, value=0.)
    AET_BSC=fipy.CellVariable(name="苔藓处的蒸发速率 (s-1)", mesh=mesh, value=0.)
    
    print("初始化边界条件...", end='->')         
    #上物理边界
    BSCfactor=BSCfactor_initial
    psi=UpBoundary.UpBoundary(mesh, X, Y, timeStepDuration,
                  psi, theta, Krel, psi_tmp, psi_max, 
                  PET, P, AET_Soil, AET_BSC,
                  PNfactor, BSCfactor,
                  psi_hygro, psi_field,                     
                  Soilfaces, BSCfaces, step=0)          
    #下物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesBottom)#隔水边界    
    #左物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesLeft)#隔水边界    
    #右物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesRight)#隔水边界
    
    # create a equation(取z轴向上为正)
    psiEqn = fipy.TransientTerm(coeff=Crel) == fipy.DiffusionTerm(coeff=Krel)+(Krel.faceValue*[[0],[1]]).divergence

    print("构建方程...", end='->') 
    #solve the equation by repeatedly looping in time:
    steps = int((Time.size)*Calibration_percent)#40%的数据用于反演参数
    desiredResidual=1e-6 
    SWC_pred=[np.mean(theta)]#初始化预测[-] 
    
    print("求解方程...", end='--->')       
    for step in range(1,steps):
        # only move forward in time once per time step
        psi.updateOld()#psi.updateOld()会自动更新所有与psi的变量        
        res = 1e+10
        MaxIterations=10
        Iterations=0
        while res > desiredResidual and Iterations < MaxIterations:            
            res = psiEqn.sweep(var=psi, dt=timeStepDuration)

            #控制约束条件
            psi.value[psi.value<psi_hygro]=psi_hygro
            psi.value[psi.value>psi_max]=psi_max
            
            #控制求解发散 
            Iterations+=1
            if Iterations>=MaxIterations: 
                print('sweep求解不收敛，解可能不准确')
                raise ValueError
            
            # print("Sweep循环中的值---------------------")
            # print("t=",step*timeStepDuration/3600,"[h]")
            # print("res=",res)
            # print("psi.value=",psi.value)
            # print("psi_tmp.value=",psi_tmp.value)
            # print("thtil=",thtil.value)
            # print("Crel.value=",Crel.value)
            # print("Krel.value=",Krel.value)         
            # print("theta=",theta.value)      
            
        print( "{:.2%}".format(step/steps), end='->')       
        #更新边界条件      
        BSCfactor=BSCfactor_final+(BSCfactor_initial-BSCfactor_final)*math.exp(-DecayRate*step)
        psi=UpBoundary.UpBoundary(mesh, X, Y, timeStepDuration,
                      psi, theta, Krel, psi_tmp, psi_max, 
                      PET, P, AET_Soil, AET_BSC,
                      PNfactor, BSCfactor, 
                      psi_hygro, psi_field,  
                      Soilfaces, BSCfaces,
                      step=step)          
        #更新上一个步长的值并保存结果       
        psi_tmp=psi.copy() #注意不能用psi_tmp=psi,psi.updateOld()会自动更新psi_tmp        
        SWC_pred.append(np.mean(theta))
        
        #图像展示每个时间步长的结果
        # theta_tmp=fipy.CellVariable(name="SWC(-)",mesh=mesh, value=theta.copy())#含水量的值
        # vi =fipy.Viewer(vars=theta_tmp, colorbar=None, datamin=0.1, datamax=0.3)
        # vi.plot()
        #fipy.TSVViewer(vars = (theta_tmp, theta_tmp.grad)).plot() 
        
    #去除数据缺失的点
    SWC_true=SWC_true[0:steps]#dataFrame数据
    SWC_true_drop=SWC_true[SWC_true.notnull()]#dataFrame数据
    SWC_true_drop=SWC_true_drop.values#dataFrame转array
    mask=np.array(SWC_pred).T#list转array 
    SWC_pred_drop=mask[~np.isnan(SWC_true.values)]
    
    #土壤含水量输出结果
    outputs=[SWC_true_drop, SWC_pred_drop, Time[0:steps], SWC_true, SWC_pred]
    # print('实测VS预测', SWC_true_drop, SWC_pred_drop,end='->')
    Error=pre_Error.Hydro_Err(true=SWC_true_drop, pred=SWC_pred_drop)
    print('预测误差', Error['rmse'], Error['nse'], Error['r_squared'], end='\n')
    return outputs

if __name__ == "__main__":
    import pandas as pd    
    
    #定义一个小函数
    def nanto1(a):
        if np.isnan(a): a=1
        return a
    
    #实验处理
    treatmentdict=['B0P0', 'B1P0', 'B2P0', 'B3P0',
                   'B0P1', 'B1P1', 'B2P1', 'B3P1', 
                   'B0P2', 'B1P2', 'B2P2', 'B3P2']
    
    #循环载入CSV文件数据
    for treatment in treatmentdict:
        file_name='{}_NSGA2_KarHydro.csv'.format(treatment)
        df=pd.read_csv(file_name)
    
        #遍历计算结果并保存成CSV
        df_result=pd.DataFrame(columns=[treatment])
        for i in list(range(-50, 0)):
            BSCfactor=nanto1(list(df['parBSCfactor'])[i])
            PNfactor=nanto1(list(df['parPNfactor'])[i])
            data=VadoseZone([BSCfactor, PNfactor], 
                            treatment=treatment, 
                            Calibration_percent=1)
            Error=pre_Error.Hydro_Err(true=data[0], pred=data[1])
            df_result['第{}次预测序列'.format(i)]=data[4]
        df_result[treatment]=data[3]
        df_result.to_csv("{}的50次预测序列".format(treatment))
    
        
        
