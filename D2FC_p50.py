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

model = Model('D2FC_p50')

### Helper Functions ###

#volume calculations
cv = 2.7e-12 #cell volume in Litres
vcm = cv/1000 #volume cell m^3
cn = 3.3 #cytoplasm to nucleus ratio
nv = cv/(cn+1) #nuc vol in L
nvm = nv/1000 #volume nucleus m^3
cytv = nv*cn #cytoplasm vol in L
pi = 3.14
r_cell = (vcm/((4*pi)/3))**(1/3) #m
r_nuc = (nvm/((4*pi)/3))**(1/3) #m
plmsa = 4*pi*(r_cell**2) #m^2
numsa = 4*pi*(r_nuc**2) #m^2
plmv = (plmsa*10e-9)*1000 # L SA in m^2 x 10 nm thickness times 1000 (m^3 to L)
numv = (numsa*10e-9)*1000 # L

AvN = 6.02214e23

Parameter('cell_vol',cv) # in Litres
Parameter('nuc_vol', nv) # in Litres
Parameter('cyt_vol', cytv) # in Litres
Parameter('cellmem_vol',plmv) # in Litres
Parameter('nucmem_vol',numv) # in Litres


#compartments
Compartment('cell_mem', dimension=2, size=cellmem_vol)
Compartment('cyt', dimension=3, size=cyt_vol, parent=cell_mem)
Compartment('nuc_mem', dimension=2, size=nucmem_vol, parent=cyt)
Compartment('nuc', dimension=3, size=nuc_vol, parent=nuc_mem)

#monomers
Monomer('IKK', ['b','s'], {'s':['n','a','i']})
Monomer('RelA',['b', 'd'])
Monomer('p50',['b', 'd'])
Monomer('IkBa', ['b','s'], {'s':['u','p']})
Monomer('tIkBa')
Monomer('A20')
Monomer('tA20')
Monomer('tTarget')
Monomer('Comp')
Monomer('tComp')

#observables
Observable('o_neut_IKKc',IKK(b=None,s='n')**cyt)
Observable('o_act_IKKc',IKK(b=None,s='a')**cyt)
Observable('o_inact_IKKc',IKK(b=None,s='i')**cyt)
Observable('o_RelA_cyt', NFkB(b=None)**cyt)
Observable('o_RelA_nuc', NFkB(b=None)**nuc)
Observable('o_IkBa_cyt', IkBa(b=None, s='u')**cyt)
Observable('o_IkBa_nuc', IkBa(b=None, s='u')**nuc)
Observable('o_IkBa_p', IkBa(b=None, s='p'))
Observable('o_IkBaNFkB_cyt', (IkBa(b=1, s='u')%NFkB(b=1))**cyt)
Observable('o_IkBaNFkB_nuc', (IkBa(b=1, s='u')%NFkB(b=1))**nuc)
Observable('o_tIkBa', tIkBa()**cyt)
Observable('o_tA20', tA20()**cyt)
Observable('o_A20', A20()**cyt)
Observable('o_tTarget', tTarget()**cyt)
Observable('o_tComp', tComp()**nuc)
Observable('o_Comp_nuc', Comp()**nuc)


#initial conditions parameters
Parameter('IkBaNFkB_0', 374353) #chose middle range value 0.3 uM from Suzan's model and converted to molecules (0.3e-6*avagodro*cytVol)
Parameter('neut_IKK_0', 100000) # .08e-6*AvN*cytv
 

#initial conditions
Initial(IKK(b=None, s='n')**cyt, neut_IKK_0)
Initial((IkBa(b=1, s='u')%NFkB(b=1))**cyt, IkBaNFkB_0)


#rate constant parameters
Parameter('kact_IKK', 0.000) # ka, set to 0.001 when TNF present, 1/s
Parameter('kinact_IKK', 0.003) #ki, 1/s
#Parameter('kneut_IKK', 0.0006) #kp, 1/s
Parameter('kbind_IkBa_NF', 8.3027e-19) # ka1a/AvN, 1/M*s*AvN IkBa + NFkB binding
Parameter('kdis_IkBa_NF', 0.05) #kd1a, 1/s IkBa:NFkB dissociation
Parameter('kphos_IkBa', 1.2288e-19) # (kc1a*1e6)/AvN, 1/M*s*AvN IKK + IkBa
Parameter('kphos_IkBaNF', 6.144e-19 ) # (kc2a*1e6)/AvN, 1/M*s*AvN IKK + IkBa:NFkB
Parameter('kdeg_pIkBa', 0.1) #kt1a and kt2a, 1/s
Parameter('kimp_NFkB', 0.0026) #ki1, 1/s
Parameter('kexp_NFkB', cn*0.000052) #ke1, 1/s
Parameter('kexp_IkBaNF', cn*0.01) #ke2a, 1/s
Parameter('kimp_IkBa', 0.00067) #ki3a, 1/s
Parameter('kexp_IkBa', cn*0.000335) #ke3a, 1/s
Parameter('ksynth_IkBa', 0.5) #c2a, 1/s
Parameter('kdeg_tIkBa', 0.0003) #c3a, 1/s
Parameter('kdeg_IkBa', 0.0005) #c4a, 1/s
Parameter('kdeg_bIkBa', 0.000022) #c5a 1/s
Parameter('ksynth_A20', 0.5) #c2 1/s
Parameter('kdeg_tA20', 0.0004) #c3 1/s
Parameter('kdeg_A20', 0.0045) #c4 1/s
Parameter('kdeg_tTarget', 0.0004) #c3 1/s
Parameter('ksynth_Comp', 0.5) #c2a 1/s
Parameter('kdeg_tComp', 0.000042857) #c6a 1/s
Parameter('kdeg_Comp', 0.0005) #c4a 1/s





#expressions
Expression('ktran_IkBa', 8.430996e10*((o_NFkB_nuc/24579)**2/((o_NFkB_nuc/24579)**2 + 1))) # Hill rate constant for IkBa transcription h=2, k1 = 0.065, kr=0, c1a = (c1a*1e-6)(AvN) 
Expression('ktran_A20', 1.204428e11*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + ((o_Comp_nuc/24579)**3) + 1))) #h =3, k1 = (0.065*1e-6)*AvN*nv, kr = 0.065, c1 = 
Expression('ktran_Target', 1.204428e11*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + (o_Comp_nuc/12289.5)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
Expression('ktran_Comp', 8.430996e10*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + (o_Comp_nuc/24579)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
 
Expression('kneut_IKK', 0.0006*(2246.12/(2246.12+o_A20)))
# =============================================================================
# Expression('ktran_IkBa', 3.702e-13*((o_NFkB_nuc/0.065)**2/((o_NFkB_nuc/0.065)**2 + 1))) # Hill rate constant for IkBa transcription h=2, k1 = 0.065, kr=0, c1a = c1a/(AvN*nv) 
# Expression('ktran_A20', 5.289e-13*((o_NFkB_nuc/0.065)**3/((o_NFkB_nuc/0.065)**3 + o_Comp_nuc*((o_Comp_nuc/0.065)**3) + 1))) #h =3, k1 = 0.065, kr = 0.065
# Expression('ktran_Target', 5.289e-13*((o_NFkB_nuc/0.065)**3/((o_NFkB_nuc/0.065)**3 + o_Comp_nuc*(o_Comp_nuc/0.065)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
# Expression('ktran_Comp', 3.702e-13*((o_NFkB_nuc/0.065)**3/((o_NFkB_nuc/0.065)**3 + o_Comp_nuc*(o_Comp_nuc/0.065)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
# 
# =============================================================================
#reaction rules
Rule('act_IKK', IKK(b=None, s='n') >> IKK(b=None, s='a'), kact_IKK) # IKK activation
Rule('inact_IKK', IKK(b=None, s='a') >> IKK(b=None, s='i'), kinact_IKK) # IKK inactivation
Rule('neut_IKK', IKK(b=None, s='i') >> IKK(b=None, s='n'), kneut_IKK) #IKK return to neutral
Rule('bind_IkBa_NFkB', IkBa(b=None, s='u') + NFkB(b=None) | 
        IkBa(b=1, s='u')%NFkB(b=1), kbind_IkBa_NF, kdis_IkBa_NF) #IkBa bind NFkB
Rule('phos_IkBa', IKK(b=None, s='a') + IkBa(b=None, s='u') >> 
     IKK(b=None, s='a') + IkBa(b=None, s='p'), kphos_IkBa) #phosphorylation IkBa
Rule('phos_IkBaNF', IKK(b=None, s='a') + IkBa(b=1, s='u')%NFkB(b=1) >> 
     IKK(b=None, s='a') + IkBa(b=1, s='p')%NFkB(b=1), kphos_IkBaNF) #phosphorylation IkBa:NFkB
Rule('deg_pIkBa', IkBa(b=None, s='p') >> None, kdeg_pIkBa) #pIkBa degradation
Rule('deg_bpIkBa', IkBa(b=1, s='p')%NFkB(b=1) >> NFkB(b=None), kdeg_pIkBa) #bound pIkBa degradation


#transport rules
Rule('import_NFkB', NFkB(b=None)**cyt | NFkB(b=None)**nuc, kimp_NFkB, kexp_NFkB) #NFkB nuclear import/export
Rule('export_IkBaNF', (IkBa(b=1, s='u')%NFkB(b=1))**nuc >> 
     (IkBa(b=1, s='u')%NFkB(b=1))**cyt, kexp_IkBaNF) #IkBa:NFkB export
Rule('import_IkBa', IkBa(b=None, s='u')**cyt | 
        IkBa(b=None, s='u')**nuc, kimp_IkBa, kexp_IkBa) #IkBa nuclear import/export

#IkBa synthesis and degradation
Rule('tran_IkBa', None >> tIkBa()**cyt, ktran_IkBa) #IkbA mRNA transcription
Rule('synth_IkBa', tIkBa()**cyt >> tIkBa()**cyt + IkBa(b=None, s='u')**cyt, ksynth_IkBa) #IkBa protein synthesis
Rule('deg_tIkBa', tIkBa() >> None, kdeg_tIkBa) #tIkBa degradation
Rule('deg_IkBa', IkBa(b=None, s='u') >> None, kdeg_IkBa) # unbound IkBa degradation
Rule('deg_bIkBa', (IkBa(b=1, s='u')%NFkB(b=1))**cyt >> NFkB(b=None), kdeg_bIkBa) # bound IkBa degradation

#A20 synthesis and degradation
Rule('tran_A20', None >> tA20()**cyt, ktran_A20) #A20 mRNA transcription
Rule('synth_A20', tA20()**cyt >> tA20()**cyt + A20()**cyt, ksynth_A20) #A20 protein synthesis
Rule('deg_tA20', tA20() >> None, kdeg_tA20) #A20 transcript degradation
Rule('deg_A20', A20() >> None, kdeg_A20) #A20 protein degradation

#prototypical inducible target
Rule('tran_Target', None >> tTarget()**cyt, ktran_Target) #inducible target transcription
Rule('deg_tTarget', tTarget() >> None, kdeg_tTarget) #tTarget degradation

#Competitor protein synthesis and degradation
Rule('tran_Comp', None >> tComp()**nuc, ktran_Comp) #Competitor mRNA transcription
Rule('synth_Comp', tComp()**nuc >> tComp()**nuc + Comp()**nuc, ksynth_Comp) #Comp protein synthesis
Rule('deg_tComp', tComp() >> None, kdeg_tComp) #Comp transcript degradation
Rule('deg_Comp', Comp() >> None, kdeg_Comp) #Comp protein degradation