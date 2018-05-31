# -*- coding: utf-8 -*-
"""
Created on Thu May 31 15:44:00 2018

@author: Amy Thurber
"""

import NFkB as m
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd

tend = 10 #time end in hourse
t = pl.linspace(0, 3600 * tend)
simres = ScipyOdeSimulator(m.model, tspan=t).run()
yout = simres.all

#set up equilibration