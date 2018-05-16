import model_Syt7at as m
import pylab as pl
import numpy as np
from pysb.integrate import *
from matplotlib import pyplot as plt

tmax = 24 #end time, in hours

t = pl.linspace(0,3600*tmax)
s = Solver(m.model,t)
s.run()

def qparams(list):
	for i in m.model.parameters:
		tag = True
		for j in list:
			if j not in i.name:
				tag = False
		if tag:
			print(i)

def chg_val(param,k):
	param.value = k
	
def graph(obs, prep_fig=True):
	s.run()
	if prep_fig:
		pl.ion()
		pl.figure()
	pl.plot(t/3600,s.yobs[obs])
	
def t_apop(apop_cutoff=5e5):
	apop_index = np.argmax(s.yobs['o_cPARP'] > apop_cutoff)
	if apop_index==0:
		if s.yobs['o_cPARP'][0] <= apop_cutoff:
			t_a = tmax
		else:
			t_a = 0
	else:
		t_a = t[apop_index]/3600.0
	return t_a				

def t_MOMP(MOMP_cutoff=1e5):
	MOMP_index = np.argmax(s.yobs['o_MOMP'] > MOMP_cutoff)
	if MOMP_index==0:
		if s.yobs['o_MOMP'][0] <= MOMP_cutoff:
			t_m = tmax
		else:
			t_m = 0
	else:
		t_m = t[MOMP_index]/3600,0
	return t_m

def create_phase():
	W = []
	X = []
	Y = []
	Z = []
	a = pl.linspace(4,8)
	A = 10**a
	for i in A:
		chg_val(m.eCa_0,i)
		s.run()
		W.append(np.log10(s.yobs['o_mCa'][-1]))
		X.append(s.yobs['o_MIMP'][-1]/1e2)
		Y.append(100*(1-s.yobs['o_ATP'][-1]/1e8))
		Z.append(s.yobs['o_iNa'][-1]/1e5)
	return W,X,Y,Z

def plot_C(param,str_range,l=0,k_str=None):
	range = []
	for i in str_range:
		range.append(float(i))
	pl.figure()
	for i in range:
		chg_val(param,i)
		W,X,Y,Z=create_phase()
		pl.plot(W,Z)
	pl.xlabel('Log Mitochondrial [Ca]')
	pl.ylabel('% Osmotic Stress')
	y = []
	if k_str is None:
		k_str=param.name
	for i in str_range:
		y.append(k_str + ' = ' + i)
	pl.legend(y,loc=l)

def plot_ATP(param,str_range,l=0,k_str=None):
    pl.figure()
    range = []
    for i in str_range:
        range.append(float(i))
    for i in str_range:
        range.append(float(i))
    for i in range:
        chg_val(param,i)
        W,X,Y,Z=create_phase()
        pl.plot(Y,Z)
    pl.xlim(0,100)
    pl.xlabel('% ATP Depletion')
    pl.ylabel('% Osmotic Stress')
    y = []
    if k_str is None:
        k_str=param.name
    for i in str_range:
        y.append(k_str + ' = ' + i)
    pl.legend(y,loc=l)

def plot_MPT(param,str_range,l=0,k_str=None):
        pl.figure()
        range = []
        for i in str_range:
                range.append(float(i))
        for i in range:
                chg_val(param,i)
                W,X,Y,Z=create_phase()
                pl.plot(X,Z)
        pl.xlabel('% MPT')
        pl.ylabel('% Osmotic Stress')
        y = []
        if k_str is None:
                k_str=param.name
        for i in str_range:
                y.append(k_str + ' = ' + i)
        pl.legend(y,loc=l)

def create_phase2(param,range):
        W = []
        X = []
        Y = []
        Z = []
        a = pl.linspace(range[0],range[1])
        A = 10**a
        for i in A:
                chg_val(param,i)
                s.run()
                X.append(s.yobs['o_MIMP'][-1]/1e2)
                Y.append(100*(1-s.yobs['o_ATP'][-1]/1e8))
                Z.append(s.yobs['o_iNa'][-1]/1e5)
        return a,X,Y,Z

def plot_A(param1,range1,param2,str_range,l=0,k_str=None):
        range = []
        for i in str_range:
                range.append(float(i))
        pl.figure()
        for i in range:
                chg_val(param2,i)
                W,X,Y,Z=create_phase2(param1,range1)
                pl.plot(W,Z)
        pl.ylabel('% Osmotic Stress')
        y = []
        if k_str is None:
                k_str=param2.name
        for i in str_range:
                y.append(k_str + ' = ' + i)
        pl.legend(y,loc=l)


def set_switch():
	chg_val(m.synth_ROS_k,0)
	W,X,Y,Z_min = create_phase2(m.PKA_0,(2,6))
	chg_val(m.synth_ROS_k,1e-1)
	W,X,Y,Z_max = create_phase2(m.PKA_0,(2,6))
	return np.array(Z_min), np.array(Z_max)

def PKA_switch(z_min, z_max, Syt7):
	chg_val(m.Syt7_0,Syt7)
	W,X,Y,Z = create_phase2(m.PKA_0,(2,6))
	z = np.array(Z)
	d_max = np.abs(z_max - z)
	d_min = np.abs(z - z_min)
	d_rel = d_max - d_min
	c_switch = W[np.argmax(d_rel>0)]
	return c_switch
