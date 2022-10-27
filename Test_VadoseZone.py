# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:07:22 2020

"""
#代码性能分析
import os
os.environ.setdefault("FIPY_SOLVERS", "scipy")
#os.environ.setdefault("FIPY_INLINE", '1')
#os.environ.setdefault("FIPY_DISPLAY_MATRIX",'terms')

import VadoseZone

import time
start =time.clock()

Parameters=[0.3, 0.5, 0.3, 0.1]
outputs=VadoseZone.VadoseZone(Parameters)
#lp.print_stats()

end = time.clock()
print('Running time: %s Seconds'%(end-start))