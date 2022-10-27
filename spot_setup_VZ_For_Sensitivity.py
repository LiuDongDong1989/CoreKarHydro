# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:01:08 2021

@author: Lenovo
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import VadoseZone_for_sensitive
import Input_Meteo
import spotpy

class spot_setup(object): 
    
    BSCfactor= Uniform(low=0.01 , high=0.99) # 0.4 #BSCfactor:苔藓覆盖蒸发系数[-]
    PNfactor=  Uniform(low=0.01 , high=0.99) #0.5  #PNfactor:松针覆盖蒸发系数[-] 
    coverage=  Uniform(low=0 , high=1) #0.4 #coverage:苔藓覆盖度[-] 
    thetas=    Uniform(low=0.2 , high=0.6)#0.60#saturation water content [-]
    thetar=    Uniform(low=0.01, high=0.1)#0.095#residual water content [-]
    alpha=     Uniform(low=0.005, high=15)#1.9#inverse of the capillary length [m-1]
    n=         Uniform(low=1.0 , high=2.0)#1.31#van Genuchten exponent [-]
    K=         Uniform(low=1.0e-08 , high=1.0e-05)#8.844444444444445e-07#Saturated hydraulic conductivity field [m/s]
    psi_initial=Uniform(low=-102/1000*3100 , high=-102/1000*30)#初始水头[m]
       
    def __init__(self, obj_func=None):
        self.obj_func = obj_func  
        
    def simulation(self, x):
        outputs = VadoseZone_for_sensitive.VadoseZone([ x[0], x[1],  x[2], x[3],  x[4], x[5], x[6], x[7], x[9] ])
        return outputs
        
    def evaluation(self):
        treatment='B1P1'
        Time, Tair, RH, P, PET, SWC_true=Input_Meteo.Meteo(0.20, 0.15, treatment=treatment)
        steps = int((Time.size))
        SWC_true=SWC_true[0:steps]
        return SWC_true
    
    def objectivefunction(self,simulation,evaluation, params=None):
        #SPOTPY expects to get one or multiple values back, 
        #that define the performance of the model run
        if not self.obj_func:
            # This is used if not overwritten by user
            like = rmse(evaluation,simulation)
        else:
            #Way to ensure flexible spot setup class
            like = self.obj_func(evaluation,simulation)    
        return like

if __name__ == "__main__":
    parallel ='seq'
    # Initialize the Hymod example
    spot_setup = spot_setup()
    
    #Select number of maximum repetitions
    # CHeck out https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/
    # How to determine an appropriate number of repetitions
    # parameters_number=9
    # rep = (1+4*16*(1+(parameters_number-2)*2))*parameters_number
    
    # #Start a sensitivity analysis
    # sampler = spotpy.algorithms.fast(spot_setup, dbname='FAST_KarHydro', dbformat='csv', db_precision=np.float32)
    # sampler.sample(rep)
    
    # Load the results gained with the fast sampler, stored in FAST_hymod.csv
    results = spotpy.analyser.load_csv_results('FAST_KarHydro')
    
    # Example plot to show the sensitivity index of each parameter
    spotpy.analyser.plot_fast_sensitivity(results, number_of_sensitiv_pars=3)
    
    # Example to get the sensitivity index of each parameter    
    SI = spotpy.analyser.get_sensitivity_of_fast(results)