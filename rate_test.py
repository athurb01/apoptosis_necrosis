# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 14:17:09 2018

@author: Amy Thurber
"""

#need to test if pysb is scaling rate parameters correctly
from pysb import *
from pysb.macros import *
import numpy as np

model = Model('rate_test')

Monomer('A',['b'])
Monomer('B',['b'])

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

#observables
Observable('o_Ac',A(b=None)**cyt)
Observable('o_Bc',B(b=None)**cyt)
Observable('o_ABc',A(b=1,)%B(b=1)**cyt)
Observable('o_An',A(b=None)**nuc)
Observable('o_Bn',B(b=None)**nuc)
Observable('o_ABn',A(b=1,)%B(b=1)**nuc)

#initial parameters
Parameter('A_0', 100000)
Parameter('B_0', 100000)

#initial conditions
Initial(A(b=None)**cyt, A_0)
Initial(B(b=None)**cyt, B_0)
Initial(A(b=None)**nuc, A_0)
Initial(B(b=None)**nuc, B_0)

#rate parameters
Parameter('kf', 8.3027e-19)
Parameter('kr', 0.05)

#rule
Rule('bind_A_B',A(b=None) + B(b=None) | A(b=1,)%B(b=1),
     kf, kr) #p65+p105 binding