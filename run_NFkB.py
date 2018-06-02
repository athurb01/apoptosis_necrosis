# -*- coding: utf-8 -*-
"""
Created on Thu May 31 15:44:00 2018

@author: Amy Thurber
"""

import NFkB as m
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd


emax = 24 #time equilibration in hours
te = pl.linspace(0, 3600 * emax)
#equilres = ScipyOdeSimulator(m.model, tspan=te).run()
#e_out = equilres.dataframe

tend = 10 #time simulation end in hours
t = pl.linspace(0, 3600 * tend)
#simres = ScipyOdeSimulator(m.model, tspan=t).run()
#y_out = simres.dataframe

#helpful functions
def chg_val(param, k):
    param.value = k

#equilibrate NFkB
def equil_NFkB(init_IKK = 100000.0, init_NFkB = 124781, run = False):
    kB_initial = [m.neut_IKK_0, m.act_IKK_0, m.inact_IKK_0, m.NFkB_0]
    kB_monomer = []
    for i in kB_initial:
        chg_val(i, 0)
        kB_monomer.append(i.name[:-2])
    
    chg_val(m.IKK_act, 0.001)
    chg_val(m.neut_IKK_0, init_IKK)
    chg_val(m.NFkB_0, init_NFkB)
    e_out = ScipyOdeSimulator(m.model, tspan=te).run().dataframe
    
    for i, initial in enumerate(kB_initial):
        chg_val(initial, e_out['o_' + kB_monomer[i]].iloc[-1])

    chg_val(m.IKK_act, 0.001)
    if run:
        simres = ScipyOdeSimulator(m.model, tspan=t).run()
        y_out = simres.dataframe

equil_NFkB()