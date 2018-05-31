# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:05:11 2018

@author: Amy Thurber
"""
# Robert has rules based model of NFkB based on the Gaudet D2FC model
# However, there seem to be some problems with how he handles compartments
# This is an attempt to create NFkB rules model using PYSB compartments 
# This model is meant to be exact match to D2FC, modeling only 1 cell using
# ka = 0.0001 when TNF present, NFkB expression = 0.3 uM = 124781 molecules

from pysb import *
from pysb.macros import *
import numpy as np

model = Model('NFkB_D2FC')

### Helper Functions ###

#Create multiple species
def create_monomers(list):
	for i in list:
		x = Monomer(i,['b'])

# E + Si <-> E:S -> E + Sa
def catalyze_act(enz,sub,site,state1,state2,kf,kr,kc):
	r1_name = '%s_%s_b' % (enz.name, sub.name)
	r2_name = '%s_%s_cat' % (enz.name, sub.name)
	E = enz(b=None)
	Si = sub({'b': None, site: state1})
	ES = enz(b=1)%sub({'b': None, site: state1})
	Sa = sub({'b': None, site: state2})
	x = Parameter(r1_name + '_kf', kf)
	y = Parameter(r1_name + '_kr', kr)
	z = Parameter(r2_name + '_k', kc)
	Rule(r1_name, E + Si | ES, x, y)
	Rule(r2_name, ES >> E + Sa, z)

# L + R <-> L:R -> Ra
def receptor_act(enz,sub,suba,kf,kr,kc):
	r1_name = '%s_%s_b' % (enz.name, sub.name)
	r2_name = '%s_%s_cat' % (enz.name, sub.name)
	E = enz(b=None)
	S = sub(b=None)
	ES = enz(b=1)%sub(b=1)
	Sa = suba(b=None)
	x = Parameter(r1_name + '_kf', kf)
	y = Parameter(r1_name + '_kr', kr)
	z = Parameter(r2_name + '_k', kc)
	Rule(r1_name, E + S | ES, x, y)
	Rule(r2_name, ES >> Sa, z)

# S -> 0
def deg(sub,k):
	r_name = 'deg_%s' % (sub.name)
	x = Parameter(r_name + '_k', k)
	Rule(r_name, sub(b=None) >> None, x)

# 0 -> S
def synth(sub,k):
	r_name = 'synth_%s' % (sub.name)
	x = Parameter(r_name + '_k', k)
	Rule(r_name, None >> sub(b=None), x)

# E -> E + S
def enz_synth(enz,prod,k):
	r_name = 'synth_%s' % (prod.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, enz() >> enz() + prod(b=None), x)

# S1 + S2 <-> S1:S2
def bind(s1,s2,kf,kr):
	r_name = '%s_%s_b' % (s1.name, s2.name)
	x = Parameter(r_name + '_kf',kf)
	y = Parameter(r_name + '_kr',kr)
	Rule(r_name, s1(b=None) + s2(b=None) | s1(b=1)%s2(b=1), x, y)

# S1 + S2 <-> S3
def bind_consume(s1,s2,s3,kf,kr):
	r_name = '%s_%s_bc' % (s1.name, s2.name)
	x = Parameter(r_name + '_kf',kf)
	y = Parameter(r_name + '_kr',kr)
	Rule(r_nam, s1(b=None) + s2(b=None) | s3(b=None), x, y)


# S1 <-> S2
def migrate_bi(s1,s2,kf,kr):
	r_name = 'mig_%s_%s' % (s1.name, s2.name)
	x = Parameter(r_name + '_kf',kf)
	y = Parameter(r_name + '_kr',kr)
	Rule(r_name, s1(b=None) | s2(b=None), x, y)

# S1 -> S2
def migrate_mono(s1,s2,k):
	r_name = 'mig_%s_%s' % (s1.name, s2.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, s1(b=None) >> s2(b=None), x)
    

#monomers
Monomer('IKK', ['b','s'], {'s':['n','a','i']})
#Monomer('p65',['b','s','l'],{'s':['u','p'],'l':['c','n']})
#Monomer('p105_p50',['s'],{'s':['u','p']})
#Monomer('p50',['b','b_DNA','b_IkB','l'],{'l':['c','n']})
#Monomer('p105',['b'])
#Monomer('NFkB',['b','s'])
#Monomer('IkB',['b','s','l'],{'s':['u','p'],'l':['c','n']})
#create_monomers(['TNF','TLR','TLRa'])

#observables
Observable('o_IKKn',IKK(b=None,s='n'))
Observable('o_IKKa',IKK(b=None,s='a'))
Observable('o_IKKi',IKK(b=None,s='i'))


#volume calculations
cv = 2.7e-12 #cell volume in Litres
vcm = cv/1000 #volume cell m^3
nc = 3.3 #cytoplasm to nucleus ratio
nv = cv/(nc+1) #nuc vol in L
nvm = nv/1000 #volume nucleus m^3
cytv = nv*nc #cytoplasm vol in L
pi = 3.14
r_cell = (vcm/((4*pi)/3))**(1/3) #m
r_nuc = (nvm/((4*pi)/3))**(1/3) #m
plmsa = 4*pi*(r_cell**2) #m^2
numsa = 4*pi*(r_nuc**2) #m^2
plmv = (plmsa*10e-9)*1000 # L SA in m^2 x 10 nm thickness times 1000 (m^3 to L)
numv = (numsa*10e-9)*1000 # L

Parameter('cell_vol',cv) # in Litres
Parameter('nuc_vol', nv) # in Litres
Parameter('cyt_vol', cytv) # in Litres
Parameter('cellmem_vol',plmv) # in Litres
Parameter('nucmem_vol',numv) # in Litres


#compartments
Compartment('cell_mem', dimension=2, size=cellmem_vol)
Compartment('cyt', dimension=3, size=cell_vol, parent=cell_mem)
Compartment('nuc_mem', dimension=2, size=nucmem_vol, parent=cyt)
Compartment('nuc', dimension=3, size=nuc_vol, parent=nuc_mem)

#initial conditions parameters
Parameter('NFkB_0', 124781) #chose middle range value 0.3 uM from Suzan's model and converted to molecules (0.3e-6*avagodro*cytVol)

#initial conditions

#rate constant parameters
Parameter('IKK_act', 0) # ka, set to 0.001 when TNF present, 1/s
Parameter('IKK_inact', 0.003) #ki, 1/s
Parameter('IKK_neut', 0.0006) #ki, 1/s

#reaction rules
Rule('act_IKK', IKK(b=None, s='n')**cyt >> IKK(b=None, s='a')**cyt, IKK_act) # IKK activation
Rule('inact_IKK', IKK(b=None, s='a')**cyt >> IKK(b=None, s='i')**cyt, IKK_inact) # IKK inactivation
Rule('inact_IKK', IKK(b=None, s='i')**cyt >> IKK(b=None, s='n')**cyt, IKK_neut) #IKK return to neutral