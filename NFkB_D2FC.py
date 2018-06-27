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

AvN = 6.02214e23

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

#monomers
Monomer('IKK', ['b','s'], {'s':['n','a','i']})
Monomer('NFkB',['b'])
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
Observable('o_NFkB_cyt', NFkB(b=None)**cyt)
Observable('o_NFkB_nuc', NFkB(b=None)**nuc)
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
Parameter('IkBaNFkB_0', 124781) #chose middle range value 0.3 uM from Suzan's model and converted to molecules (0.3e-6*avagodro*cytVol)
Parameter('neut_IKK_0', 100000) # .08e-6*AvN*cytv
Parameter('act_IKK_0', 0)
Parameter('inact_IKK_0', 0) 

#initial conditions
Initial(IKK(b=None, s='n')**cyt, neut_IKK_0)
Initial((IkBa(b=1, s='u')%NFkB(b=1))**cyt, IkBaNFkB_0)


#rate constant parameters
Parameter('kact_IKK', 0.000) # ka, set to 0.001 when TNF present, 1/s
Parameter('kinact_IKK', 0.003) #ki, 1/s
Parameter('kneut_IKK', 0.0006) #ki, 1/s
Parameter('kbind_IkBa_NF', 8.3027e-19) # ka1a/AvN, 1/M*s*AvN IkBa + NFkB binding
Parameter('kdis_IkBA_NF', 0.05) #kd1a, 1/s IkBa:NFkB dissociation
Parameter('kphos_IkBa', 1.2288e-19) # (kc1a*1e6)/AvN, 1/M*s*AvN IKK + IkBa
Parameter('kphos_IkBaNF', 6.144e-19 ) # (kc2a*1e6)/AvN, 1/M*s*AvN IKK + IkBa:NFkB
Parameter('kdeg_pIkBa', 0.1) #kt1a, 1/s
Parameter('kimp_NFkB', 0.0026) #ki1, 1/s
Parameter('kexp_NFkB', 0.000052) #ki1, 1/s
Parameter('kexp_IkBaNF', 0.01) #ke2a, 1/s
Parameter('kimp_IkBa', 0.00067) #ki3a, 1/s
Parameter('kexp_IkBa', 0.000335) #ke3a, 1/s
Parameter('ksynth_IkBa', 0.5) #c2a, 1/s
Parameter('kdeg_tIkBa', 0.0003) #c3a, 1/s
Parameter('kdeg_IkBa', 0.0005) #c4a, 1/s
Parameter('kdeg_bIkBa', 0.000022) #c5a 1/s
Parameter('ksynth_A20', 0.5) #c2 1/s
Parameter('kdeg_tA20', 0.0004) #c3 1/s
Parameter('kdeg_A20', 0.0045) #c4 1/s
Parameter('kdeg_tTarget', 0.0004) #c3 1/s
Parameter('ksynth_Comp', 0.5) #c2a 1/s
Parameter('kdeg_tComp', 0.00004) #c6a 1/s
Parameter('kdeg_Comp', 0.0005) #c4a 1/s

#expressions
Expression('ktran_IkBa', 0.05294*((o_NFkB_nuc/24579)**2/((o_NFkB_nuc/24579)**2 + 1))) # Hill rate constant for IkBa transcription h=2, k1 = 0.065, kr=0, c1a = c1a/(AvN*nv) 
Expression('ktran_A20', 0.07563*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + o_Comp_nuc*((o_Comp_nuc/24579)**3) + 1))) #h =3, k1 = 0.065, kr = 0.065
Expression('ktran_Target', 0.07563*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + o_Comp_nuc*(o_Comp_nuc/24579)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
Expression('ktran_Comp', 0.05294*((o_NFkB_nuc/24579)**3/((o_NFkB_nuc/24579)**3 + o_Comp_nuc*(o_Comp_nuc/24579)**3 + 1))) #h =3, k1 = 0.065, kr = 0.065
 
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
        IkBa(b=1, s='u')%NFkB(b=1), kbind_IkBa_NF, kdis_IkBA_NF) #IkBa bind NFkB
Rule('phos_IkBa', IKK(b=None, s='a') + IkBa(b=None, s='u') >> 
     IKK(b=None, s='a') + IkBa(b=None, s='p'), kphos_IkBa) #phosphorylation IkBa
Rule('phos_IkBaNF', IKK(b=None, s='a') + IkBa(b=1, s='u')%NFkB(b=1) >> 
     IKK(b=None, s='a') + IkBa(b=1, s='p')%NFkB(b=1), kphos_IkBaNF) #phosphorylation IkBa:NFkB
Rule('deg_pIkBa', IkBa(s='p') >> None, kdeg_pIkBa) #pIkBa degradation

#transport rules
Rule('import_NFkB', NFkB(b=None)**cyt | NFkB(b=None)**nuc, kimp_NFkB, kexp_NFkB) #NFkB nuclear import/export
Rule('export_IkBaNF', (IkBa(b=1, s='u')%NFkB(b=1))**nuc >> 
     (IkBa(b=1, s='u')%NFkB(b=1))**cyt, kexp_IkBaNF) #IkBa:NFkB export
Rule('import_IkBa', IkBa(b=None, s='u')**cyt | 
        IkBa(b=None, s='u')**nuc, kimp_IkBa, kexp_IkBa) #IkBa nuclear import/export

#IkBa synthesis and degradation
Rule('tran_IkBa', NFkB(b=None)**nuc >> NFkB(b=None)**nuc + tIkBa()**cyt, ktran_IkBa) #IkbA mRNA transcription
Rule('synth_IkBa', tIkBa()**cyt >> tIkBa()**cyt + IkBa(b=None, s='u')**cyt, ksynth_IkBa) #IkBa protein synthesis
Rule('deg_tIkBa', tIkBa() >> None, kdeg_tIkBa) #tIkBa degradation
Rule('deg_IkBa', IkBa(b=None, s='u')**cyt >> None, kdeg_IkBa) # unbound IkBa degradation
Rule('deg_bIkBa', (IkBa(b=1, s='u')%NFkB(b=1))**cyt >> NFkB(b=None), kdeg_bIkBa) # bound IkBa degradation

#A20 synthesis and degradation
Rule('tran_A20', NFkB(b=None)**nuc >> NFkB(b=None)**nuc + tA20()**cyt, ktran_A20) #A20 mRNA transcription
Rule('synth_A20', tA20()**cyt >> tA20()**cyt + A20()**cyt, ksynth_A20) #A20 protein synthesis
Rule('deg_tA20', tA20() >> None, kdeg_tA20) #A20 transcript degradation
Rule('deg_A20', A20() >> None, kdeg_A20) #A20 protein degradation

#prototypical inducible target
Rule('tran_Target', NFkB(b=None)**nuc >> NFkB(b=None)**nuc + tTarget()**cyt, ktran_Target) #inducible target transcription
Rule('deg_tTarget', tTarget() >> None, kdeg_tTarget) #tTarget degradation

#Competitor protein synthesis and degradation
Rule('tran_Comp', NFkB(b=None)**nuc >> NFkB(b=None)**nuc + tComp()**nuc, ktran_Comp) #Competitor mRNA transcription
Rule('synth_Comp', tComp()**nuc >> tComp()**nuc + Comp()**nuc, ksynth_Comp) #Comp protein synthesis
Rule('deg_tComp', tComp() >> None, kdeg_Comp) #Comp transcript degradation
Rule('deg_Comp', Comp() >> None, kdeg_Comp) #Comp protein degradation