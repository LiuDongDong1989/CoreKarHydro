# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 16:21:12 2021

@author: Lenovo

This example implements the python version of hymod into SPOTPY.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import Input_Meteo
import spotpy
import matplotlib.pyplot as plt
import numpy as np
        
def evaluation(treatment):   
    Time, Tair, RH, P, PET, SWC_true=Input_Meteo.Meteo(0.20, 0.15, treatment=treatment)
    steps = int((Time.size)*0.40)#40%的数据用于反演参数
    SWC_true=SWC_true[0:steps]#dataFrame数据
    SWC_true_drop=SWC_true[SWC_true.notnull()]#dataFrame数据
    SWC_true_drop=SWC_true_drop.values#dataFrame转array  
    return SWC_true_drop
        
if __name__ == "__main__":
    
    treatment='B3P2'
    csv_name='{}_NSGA2_KarHydro'.format(treatment)    
    results = spotpy.analyser.load_csv_results(csv_name) 
    fontsize_xyl=24
    
    #画目标函数的变化
    fig= plt.figure(1,figsize=(16,9),dpi=300)
    plt.plot(results['like1'])
    plt.title('{}'.format(treatment),fontsize=fontsize_xyl)
    plt.ylabel('RMSE', fontsize=fontsize_xyl)
    plt.xlabel('Iteration', fontsize=fontsize_xyl)
    plt.xlim(0, 5000)
    plt.ylim()
    plt.xticks(fontsize=fontsize_xyl)
    plt.yticks(fontsize=fontsize_xyl)
    fig.savefig('{}_objectivefunction.png'.format(csv_name), dpi=600)
    
    #画反演的不确定性#
    fields=[word for word in results.dtype.names if word.startswith('sim')]    
    q5,q25,q75,q95=[],[],[],[]
    for field in fields:
        q5.append(np.percentile(results[field][-100:-1],2.5))
        q95.append(np.percentile(results[field][-100:-1],97.5))
        
    fig= plt.figure(2, figsize=(16,9), dpi=300)    
    plt.plot(q5,color='dimgrey',linestyle='solid')
    plt.plot(q95,color='dimgrey',linestyle='solid')
    plt.fill_between(np.arange(0,len(q5),1),list(q5),list(q95),
                      facecolor='dimgrey',zorder=0,
                      linewidth=0, label='Simulated')  
    plt.plot(evaluation(treatment), 'ro', label='Observed', markersize=16)
    
    plt.title('{}'.format(treatment),fontsize=fontsize_xyl)
    plt.xlabel("Evaporation Days", fontsize=fontsize_xyl)
    plt.ylabel("SWC (-)", fontsize=fontsize_xyl)
    plt.xlim(0, 25)
    plt.ylim(0.25,0.32)
    plt.xticks(np.arange(0,    25,   step=2),    fontsize=fontsize_xyl)
    plt.yticks(np.arange(0.25, 0.35, step=0.01), fontsize=fontsize_xyl)
    plt.legend(fontsize=fontsize_xyl, loc='upper right')
    
    fig.savefig('{}_Simulated.png'.format(csv_name), dpi=600)
    
    