# -*- coding: utf-8 -*-
"""
Created on Thu May 31 15:44:00 2018

@author: Amy Thurber
"""

import D2FC_p50 as m
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd
import numpy as np
## to re-import an updated model m.model.reload()
AvN = 6.02214e23
emax = 24 #time equilibration in hours
te = pl.linspace(0, 3600 * emax)
tend = 10 #time simulation end in hours
t = pl.linspace(0, 3600 * tend)

#helpful dictionaries of parameters
equil_dict = {'kact_IKK':0}
steady_protein_dict = {'kact_IKK':0, 'ksynth_IkBa':0, 'kdeg_IkBa':0,
                       'kdeg_bIkBa':0, 'ksynth_A20':0, 'kdeg_A20':0,
                       'ksynth_p50':0, 'kdeg_p50':0}

#helpful lists of observables
monomers = ['o_RelAm_cyt','o_RelAm_nuc', 'o_p50m_cyt', 'o_p50m_nuc']
dimers_free = ['o_A50_cyt', 'o_A50_nuc','o_AA_cyt', 'o_AA_nuc', 'o_p5050_cyt', 'o_p5050_nuc' ]
IkBa = ['o_IkBa_cyt', 'o_IkBa_nuc', 'o_IkBa_p', 'o_IkBab_p', 'o_IkBaA50_cyt',
        'o_IkBaA50_nuc', 'o_IkBap5050_cyt', 'o_IkBap5050_nuc']
p50 = []
transcripts = ['o_tIkBa', 'o_tp50', 'o_tTarget',  'o_tA20']
# to equilibrate and use equilibration values as initial for simulation
sim = ScipyOdeSimulator(m.model, te)
eq_res = sim.run(param_values = steady_protein_dict)
df_equil = eq_res.dataframe
#add addition columns that require calculations
df_equil['o_ratio_nuc'], df_equil['o_ratio_cyt'] = [
        df_equil['o_A50_nuc']/df_equil['o_p5050_nuc'],
        df_equil['o_A50_cyt']/df_equil['o_p5050_cyt']]
df_equil['o_p50_total'] = (df_equil['o_p50m_nuc'] + df_equil['o_p50m_cyt'] + 
    df_equil['o_A50_nuc'] + df_equil['o_A50_cyt'] + 
    2*(df_equil['o_p5050_nuc'] + df_equil['o_p5050_cyt']) + 
    df_equil['o_IkBaA50_nuc'] + df_equil['o_IkBaA50_cyt'] + 
    2*(df_equil['o_IkBap5050_nuc'] + df_equil['o_IkBap5050_cyt']))
df_equil['o_RelA_total'] = (df_equil['o_RelAm_nuc'] + df_equil['o_RelAm_cyt'] + 
    df_equil['o_A50_nuc'] + df_equil['o_A50_cyt'] + 
    2*(df_equil['o_AA_nuc'] + df_equil['o_AA_cyt']) + 
    df_equil['o_IkBaA50_nuc'] + df_equil['o_IkBaA50_cyt'] + 
    2*(df_equil['o_IkBaAA_nuc'] + df_equil['o_IkBaAA_cyt']))

df_equil.plot() #will plot all species
eq = df_equil.iloc[-1, 0:len(m.model.species)] #grabs species and not obervables

act_dict = {'kact_IKK': 0.001}
sim_res = sim.run(initials = eq, param_values = param_dict, tspan = t)
df_sim = sim_res.dataframe
df_sim.plot()

## to reset index so sim time continues from end equilibration
#df.index + emax
# pd.concat(df_equil, df_sim)

##things I'm not ready to delete
#equilres = ScipyOdeSimulator(m.model, tspan=te).run()
#e_out = equilres.dataframe
# sim = ScipyOdeSimulator(m.model, t)
#y_out = sim.run().dataframe

#parameter sweep function
def param_sweep(param, low, high, outputs):
    sim = ScipyOdeSimulator(m.model, te)
    results = []
    params = {'k1': 0.03, 'k2': 10000}
    for k in np.linspace(low, high, 10):
        params['k3'] = k3
        result = sim.run(param_values=params)
        results.append(result)

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