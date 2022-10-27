# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 13:01:12 2019
---------------------函数说明-------------------------------------------------
计算预测误差
-----------------------------------------------------------------------------
"""
from sklearn.metrics import max_error, mean_absolute_error, mean_squared_error, mean_squared_log_error, median_absolute_error, r2_score
import HydroErr as he
import numpy as np

def sklearn_Err(true=0, pred=0):
    pre_Error={}
    pre_Error['MaxError'] = max_error(true, pred)
    pre_Error['MAE'] = mean_absolute_error(true, pred)
    pre_Error['MSE'] = mean_squared_error(true, pred)
    pre_Error['MSLE'] = mean_squared_log_error(true, pred)
    pre_Error['MedAE'] = median_absolute_error(true, pred)
    pre_Error['R2'] = r2_score(true, pred)
    return pre_Error

def Hydro_Err(true=0, pred=0):
    sim=pred
    obs=true
    pre_Error={}#目标函数
    optimalValue={}#各目标函数最优值
    worstValue={}#各目标函数最差值
    distanceValue={}#距离函数最差值

    pre_Error['me']=he.me(sim, obs)
    optimalValue['me']=0
    worstValue['me']=[-float('inf'), float('inf')]

    pre_Error['mae']=he.mae(sim, obs)
    optimalValue['mae']=0
    worstValue['mae'] =float('inf')
    
    pre_Error['mse']=he.mse(sim, obs)
    optimalValue['mse']=0
    worstValue['mse'] = float('inf')
    
    pre_Error['mle']=he.mse(sim, obs)
    optimalValue['mle']=0
    worstValue['mle'] = [-float('inf'),float('inf')]
    
    pre_Error['male']=he.male(sim, obs)
    optimalValue['male']=0
    worstValue['male'] = float('inf')
    
    pre_Error['msle']=he.msle(sim, obs)
    optimalValue['msle']=0
    worstValue['msle'] = float('inf')
    
    pre_Error['mde']=he.mde(sim, obs)
    optimalValue['mde']=0
    worstValue['mde'] = [-float('inf'),float('inf')]
    
    pre_Error['mdae']=he.mdae(sim, obs)
    optimalValue['mdae'] = 0
    worstValue['made'] = float('inf')

    pre_Error['mdse'] = he.mdse(sim, obs)
    optimalValue['mdse'] = 0
    worstValue['mdse'] = float('inf')

    pre_Error['ed'] = he.ed(sim, obs)
    optimalValue['ed'] = 0
    worstValue['ed'] = float('inf')

    pre_Error['ned'] = he.ned(sim, obs)
    optimalValue['ned'] = 0
    worstValue['ned'] = float('inf')

    pre_Error['rmse'] = he.rmse(sim, obs)
    optimalValue['rmse'] = 0
    worstValue['rmse'] = float('inf')

    pre_Error['rmsle'] = he.rmsle(sim, obs)
    optimalValue['rmsle'] = 0
    worstValue['rmsle'] = float('inf')

    pre_Error['nrmse_range']=he.nrmse_range(sim, obs)
    optimalValue['nrmse_range']=0
    worstValue['nrmse_range'] = float('inf')

    pre_Error['nrmse_mean'] = he.nrmse_mean(sim, obs)
    optimalValue['nrmse_mean'] = 0
    worstValue['nrmse_mean'] = float('inf')

    pre_Error['nrmse_iqr'] = he.nrmse_iqr(sim, obs)
    optimalValue['nrmse_iqr'] = 0
    worstValue['nrmse_iqr'] = float('inf')

    pre_Error['irmse'] = he.irmse(sim, obs)
    optimalValue['irmse'] = 0
    worstValue['irmse'] = float('inf')

    pre_Error['mase'] = he.mase(sim, obs)
    optimalValue['mase'] =0
    worstValue['mase'] =float('inf')

    pre_Error['r_squared'] = he.r_squared(sim, obs)
    optimalValue['r_squared'] = 1
    worstValue['r_squared'] = 0

    pre_Error['pearson_r'] = he.pearson_r(sim, obs)
    optimalValue['pearson_r'] = 1
    worstValue['pearson_r'] = -1

    pre_Error['spearman_r'] = he.spearman_r(sim, obs)
    optimalValue['spearman_r'] = 1
    worstValue['spearman_r'] = -1

    pre_Error['acc'] = he.acc(sim, obs)
    optimalValue['acc'] = 1
    worstValue['acc'] = -1

    pre_Error['mape'] = he.mape(sim, obs)
    optimalValue['mape'] = 0
    worstValue['mape'] =float('inf')

    pre_Error['mapd']=he.mapd(sim, obs)
    optimalValue['mapd']=0
    worstValue['mapd'] =float('inf')

    pre_Error['maape'] = he.maape(sim, obs)
    optimalValue['maape'] = 0
    worstValue['maape'] = np.pi/2

    pre_Error['smape1'] = he.smape1(sim, obs)
    optimalValue['smape1'] = 0
    worstValue['smape1'] = 1

    pre_Error['smape2'] = he.smape2(sim, obs)
    optimalValue['smape2'] = 0
    worstValue['smape2'] = 2

    pre_Error['d'] = he.d(sim, obs)
    optimalValue['d'] = 1
    worstValue['d'] = 0

    pre_Error['d1'] = he.d1(sim, obs)
    optimalValue['d1'] = 1
    worstValue['d1'] = 0

    pre_Error['dmod'] = he.dmod(sim, obs)
    optimalValue['dmod'] = 1
    worstValue['dmod'] = 0

    pre_Error['drel'] = he.drel(sim, obs)
    optimalValue['drel'] = 1
    worstValue['drel'] = 0

    pre_Error['dr'] = he.dr(sim, obs)
    optimalValue['dr'] = 1
    worstValue['dr'] = -1

    pre_Error['watt_m'] = he.watt_m(sim, obs)
    optimalValue['watt_m'] = 1
    worstValue['watt_m'] = -1

    pre_Error['mb_r'] = he.mb_r(sim, obs)
    optimalValue['mb_r'] = 1
    worstValue['mb_r'] = 0

    pre_Error['nse'] = he.nse(sim, obs)
    optimalValue['nse'] = 1
    worstValue['nse'] = -float('inf')

    pre_Error['nse_mod'] = he.nse_mod(sim, obs)
    optimalValue['nse_mod'] = 1
    worstValue['nse_mod'] = -float('inf')

    pre_Error['nse_rel'] = he.nse_rel(sim, obs)
    optimalValue['nse_rel'] = 1
    worstValue['nse_rel'] = -float('inf')

    pre_Error['kge_2009'] = he.kge_2009(sim, obs)
    optimalValue['kge_2009'] = 1
    worstValue['kge_2009'] = -float('inf')

    pre_Error['kge_2012'] = he.kge_2012(sim, obs)
    optimalValue['kge_2012'] = 1
    worstValue['kge_2012'] = -float('inf')

    pre_Error['lm_index'] = he.lm_index(sim, obs)
    optimalValue['lm_index'] = 1
    worstValue['lm_index'] = 0

    pre_Error['d1_p'] = he.d1_p(sim, obs)
    optimalValue['d1_p'] = 1
    worstValue['d1_p'] = 0

    pre_Error['ve'] = he.ve(sim, obs)
    optimalValue['ve'] = 0
    worstValue['ve'] = 1

    pre_Error['sa'] = he.sa(sim, obs)
    optimalValue['sa'] = 0
    worstValue['sa'] = [-np.pi/2,np.pi/2]

    pre_Error['sc'] = he.sc(sim, obs)
    optimalValue['sc'] = 0
    worstValue['sc'] = [-np.pi/2,np.pi/2]

    pre_Error['sid'] = he.sid(sim, obs)
    optimalValue['sid'] = 0
    worstValue['sid'] = [-np.pi/2,np.pi/2]

    pre_Error['sga'] = he.sga(sim, obs)
    optimalValue['sga'] = 0
    worstValue['sga'] = [-np.pi/2,np.pi/2]
   
    for key in pre_Error:
        if pre_Error[key] >= optimalValue [key]:
            distanceValue[key] = pre_Error[key] -optimalValue[key]
        else:
            distanceValue[key] = optimalValue[key] - pre_Error[key]
            
    # print(pre_Error)      
    #返回三个字典[实际目标函数，目标函数的最优值，目标函数的最差值, 离最优值的距离]
    return pre_Error
    
if __name__ == "__main__":
    pred=[0.24154472, 0.15332635, 0.19890591, 0.28256742, 0.21366198,0.49954502, 0.21811334, 0.14823118]
    true=[0.24154472, 0.19037044, 0.24516848, 0.26783164, 0.20595483,0.36159936, 0.19061863, 0.25961065]
    Error1=Hydro_Err(true, pred)
    Error2=sklearn_Err(true,pred)
    print(Error1)
    print(Error2)
    
    
    
    
    
    
    
