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

import os
os.environ.setdefault("FIPY_SOLVERS", "scipy")

from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import VadoseZone
import Input_Meteo
import spotpy
import matplotlib.pyplot as plt
import numpy as np

class spot_setup(object):
    treatment='B0P1'
    if treatment[1]=='0' and treatment[-1]=='0':
        thetas=Uniform(low=0.3, high=0.6)
        thetar=Uniform(low=0.0, high=0.2)
        alpha=Uniform(low=1, high=6)
        n=Uniform(low=1.1, high=1.7)
        Ks=Uniform(low=1e-8, high=1e-4)
    elif treatment[1]!='0' and treatment[-1]!='0':
        BSCfactor = Uniform(low=0.3, high=1)
        PNfactor  = Uniform(low=0.3, high=1)
        DecayFactor=Uniform(low=0.3, high=1)
        DecayRate = Uniform(low=0.0, high=0.5)
    elif treatment[1]!='0' and treatment[-1]=='0':
        BSCfactor = Uniform(low=0.3, high=1)
        DecayFactor=Uniform(low=0.3, high=1)
        DecayRate = Uniform(low=0.0, high=0.5)
    elif treatment[1]=='0' and treatment[-1]!='0':
        PNfactor  = Uniform(low=0.3, high=1)
      
    def __init__(self, treatment=treatment, obj_func=None):
        self.obj_func = obj_func  
        self.treatment=treatment
        
    def simulation(self, x):
        if self.treatment[1]=='0' and self.treatment[-1]=='0':
            outputs = VadoseZone.VadoseZone([x[0], x[1], x[2], x[3], x[4] ], self.treatment)  
        elif self.treatment[1]!='0' and self.treatment[-1]!='0':
            outputs = VadoseZone.VadoseZone([ x[0], x[1], x[2], x[3] ], self.treatment) 
        elif self.treatment[1]!='0' and self.treatment[-1]=='0':
            outputs = VadoseZone.VadoseZone([x[0], x[1], x[2]], self.treatment)
        elif self.treatment[1]=='0' and self.treatment[-1]!='0':
            outputs = VadoseZone.VadoseZone([x[0]], self.treatment)
        return outputs[1]
        
    def evaluation(self):   
        Time, Tair, RH, P, PET, SWC_true=Input_Meteo.Meteo(0.20, 0.15, treatment=self.treatment)
        steps = int((Time.size)*0.40)#40%的数据用于反演参数
        SWC_true=SWC_true[0:steps]#dataFrame数据
        SWC_true_drop=SWC_true[SWC_true.notnull()]#dataFrame数据
        SWC_true_drop=SWC_true_drop.values#dataFrame转array  
        return SWC_true_drop
    
    def objectivefunction(self, simulation, evaluation, params=None):
        #SPOTPY expects to get one or multiple values back, 
        #that define the performance of the model run
        if not self.obj_func:
            # This is used if not overwritten by user
            like = rmse(evaluation,simulation)
        else:
            #Way to ensure flexible spot setup class
            like = self.obj_func(evaluation,simulation)    
        return like

def multi_obj_func(evaluation, simulation):
    #used to overwrite objective function
    like1 = spotpy.objectivefunctions.rmse(evaluation, simulation)
    like2 = spotpy.objectivefunctions.nashsutcliffe(evaluation, simulation)*-1
    like3 = spotpy.objectivefunctions.rsquared(evaluation, simulation)*-1
    print("目标函数是", [like1, like2, like3])
    return np.array([like1, like2])
    
if __name__ == "__main__":

    generations=5
    n_pop = 10
    skip_duplicates = True
    
    sp_setup=spot_setup( obj_func=multi_obj_func)
    sampler = spotpy.algorithms.NSGAII(spot_setup=sp_setup, dbname='{}_NSGA2_KarHydro'.format(sp_setup.treatment), dbformat="csv") 
    
    sampler.sample(generations, n_obj= 2, n_pop = n_pop)
    
    results = spotpy.analyser.load_csv_results('{}_NSGA2_KarHydro'.format(sp_setup.treatment))
       
    fig= plt.figure(1,figsize=(9,5))
    plt.plot(results['like2'])
    plt.show()
    plt.ylabel('Objective function')
    plt.xlabel('Iteration')
    fig.savefig('{}_Objectivefunction.png'.format(sp_setup.treatment),dpi=300)
    
    # Example plot to show remaining parameter uncertainty #
    fields=[word for word in results.dtype.names if word.startswith('sim')]
    fig= plt.figure(3, figsize=(16,9))
    ax11 = plt.subplot(1,1,1)
    q5,q25,q75,q95=[],[],[],[]
    for field in fields:
        q5.append(np.percentile(results[field][-100:-1],2.5))
        q95.append(np.percentile(results[field][-100:-1],97.5))
    ax11.plot(q5,color='dimgrey',linestyle='solid')
    ax11.plot(q95,color='dimgrey',linestyle='solid')
    ax11.fill_between(np.arange(0,len(q5),1),list(q5),list(q95),facecolor='dimgrey',zorder=0,
                    linewidth=0,label='parameter uncertainty')  
    ax11.plot(sp_setup.evaluation(),'r.',label='data')
    ax11.legend()
    fig.savefig('{}_NSGA2_KarHydro.png'.format(sp_setup.treatment),dpi=300)
    #########################################################