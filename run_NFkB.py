# -*- coding: utf-8 -*-
"""
Created on Thu May 31 15:44:00 2018

@author: Amy Thurber
"""

import NFkB_D2FC as m
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd
## to re-import an updated model model.reload()

emax = 24 #time equilibration in hours
te = pl.linspace(0, 3600 * emax)
#equilres = ScipyOdeSimulator(m.model, tspan=te).run()
#e_out = equilres.dataframe

tend = 10 #time simulation end in hours
t = pl.linspace(0, 3600 * tend)
# sim = ScipyOdeSimulator(m.model, t)
#y_out = sim.run().dataframe

# to equilibrate and use equilibration values as initial for simulation
sim = ScipyOdeSimulator(m.model, te)
# no IkBa degradation
equil_dict = {'kact_IKK':0}
eq_res = sim.run(param_values = equil_dict)
df_equil = eq_res.dataframe
df_equil.plot() #will plot all species
eq = df_equil.iloc[-1, 0:len(m.model.species)] #grabs species and not obervables
param_dict = {'kact_IKK': 0.001}
sim_res = sim.run(initials = eq, param_values = param_dict, tspan = t)
df_sim = sim_res.dataframe
df_sim.plot()

## to reset index so sim time continues from end equilibration
#df.index + emax
# pd.concat(df_equil, df_sim)


#helpful functions
def chg_val(param, k):
    param.value = k

#equilibrate NFkB
def equil_NFkB(init_IKK = 100000.0, init_NFkB = 124781, sim_vals, run = False):
    equil_vals = {'IKK_act':0.001, 'neut_IKK_0':init_IKK, 'act_IKK_0':0,
                       'inact_IKK_0':0, 'NFkB_0':init_NFkB}
    e_out = sim.run(param_values = equil_vals).dataframe
    for k, value in equil_vals:
        init_param = m.
    
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