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

def VadoseZone(Parameters):
    
    BSCfactor=Parameters[0] #0.4 #BSCfactor:苔藓覆盖蒸发系数[-]
    PNfactor=Parameters[1]  #0.5  #PNfactor:松针覆盖蒸发系数[-] 
    coverage=Parameters[2]  #0.4 #coverage:苔藓覆盖度[-] 
    thetas=Parameters[3]    #0.60#saturation water content [-]
    thetar=Parameters[4]    #0.095#residual water content [-]
    alpha=Parameters[5]     #1.9#inverse of the capillary length [m-1]
    n=Parameters[6]         #1.31#van Genuchten exponent [-]
    K=Parameters[7]         #8.844444444444445e-07#Saturated hydraulic conductivity field [m/s]
    psi_initial=Parameters[8]#初始水头[m]
    
    print('苔藓覆盖蒸发系数=',BSCfactor)
    print('松针覆盖蒸发系数=',PNfactor)
          
    import fipy.solvers.solver 
    print("判断求解器是否采用petsc求解器", end='/*/')
    print('调用求解器...', fipy.solvers.solver)
      
    print("划分网格...", end='/*/')
    nx, ny =10, 10
    Lx, Ly = 0.20,0.15
    dx, dy = Lx/nx, Ly/ny #单位：m
    mesh = fipy.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
    timeStepDuration=86400#[s]#timeStepDuration =0.9* dy**2 / (2 * max(Krel))
    X, Y = mesh.faceCenters#单元面的坐标
    
    print("定义覆盖区域、输入气象数据...", end='/*/')
    Time, Tair, RH, P, PET, SWC_true=Input_Meteo.Meteo(Lx, Ly)
    Soilfaces, BSCfaces=Define_Cover.Define_Cover(mesh, X, Y, Lx, Ly, coverage=coverage)
    
    print("初始化水头分布...", end='/*/')
    # Water pressure field [m] - field of resolution
    psi_max=0.00 #地表最大滞蓄深度[m]
    #psi_initial=-10  #初始水头[m]
    psi=fipy.CellVariable(name="Total water pressure head (m)", mesh=mesh,value=psi_initial, hasOld=1)#当前时间步长的值，扩散系数是非线性的,用hasOld来保持原值直至更新动作
    psi_tmp = fipy.CellVariable(name="Total water pressure head of last time step (m)", mesh=mesh, value=psi_initial)#上一个时间步长的值 
    
    print("定义土壤水力特征参数...", end='/*/')
    psi_hygro=-102/1000*3100 #土壤吸湿点[m]  
    psi_field=-102/1000*30#土壤田持含水量下的水吸力[m] 
    thtil, Krel, Crel, theta=SWRC.VGmodel(psi, thetas=thetas, thetar=thetar, alpha=alpha, n=n, S=0.00001, K= K)

    print("计算潜水蒸发、初始化蒸发...", end='/*/')
    # from Calculation_PET import Calculation_PET
    # PET_pre=Calculation_PET(Tair=Tair, RH=RH, Lx=Lx)#大气蒸发量的预测值
    AET_Soil=fipy.CellVariable(name="裸土处的蒸发速率 (s-1)", mesh=mesh, value=0.)
    AET_BSC=fipy.CellVariable(name="苔藓处的蒸发速率 (s-1)", mesh=mesh, value=0.)
    
    print("初始化边界条件...", end='/*/')        
    #上物理边界
    psi=UpBoundary.UpBoundary(mesh, X, Y, timeStepDuration,
                  psi, theta, Krel, psi_tmp, psi_max, 
                  PET, P, AET_Soil, AET_BSC,
                  PNfactor, BSCfactor,
                  psi_hygro, psi_field,                     
                  Soilfaces, BSCfaces,
                  step=0)               
    #下物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesBottom)#隔水边界    
    #左物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesLeft)#隔水边界    
    #右物理边界
    psi.faceGrad.constrain([[0],[0]], mesh.facesRight)#隔水边界
     
    print("构建方程...", end='/*/')    
    # create a equation(取z轴向上为正)
    psiEqn = fipy.TransientTerm(coeff=Crel) == fipy.DiffusionTerm(coeff=Krel)+(Krel.faceValue*[[0],[1]]).divergence    
    steps = int((Time.size))
    desiredResidual=1e-6 
    SWC_pred=[np.mean(theta)]#初始化预测[-] 
    
    print("求解方程...", end='/*/')   
    for step in range(1,steps):#solve the equation by repeatedly looping in time
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
        
        print( "{:.2%}".format(step/steps), end='--->' )  
        #更新边界条件       
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
    
    print('保存结果')        
    #土壤含水量输出结果
    outputs=SWC_pred
      
    return outputs

if __name__ == "__main__":
    import time
    start =time.clock()
    VadoseZone([0.4, 0.5,0.4, 0.6, 0.095, 1.9, 1.31, 8.844444444444445e-07])    
    end = time.clock()
    print('Running time: %s Seconds'%(end-start))
