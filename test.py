# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 13:47:10 2018

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

model = Model('test')

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
Monomer('RelA', ['b', 'd'])
Monomer('p50', ['b', 'd'])
Monomer('IkBa', ['b','s'], {'s':['u','p']})
Monomer('tIkBa')
Monomer('A20')
Monomer('tA20')
Monomer('tTarget')
Monomer('tp50')

#observables
Observable('o_neut_IKKc', IKK(b=None,s='n')**cyt)
Observable('o_act_IKKc', IKK(b=None,s='a')**cyt)
Observable('o_inact_IKKc',IKK(b=None,s='i')**cyt)
Observable('o_RelAm_cyt', RelA(b=None, d=None)**cyt)
Observable('o_RelAm_nuc', RelA(b=None, d=None)**nuc)
Observable('o_p50m_nuc', p50(b=None, d=None)**nuc)
Observable('o_p50m_cyt', p50(b=None, d=None)**cyt)
Observable('o_A50_cyt', RelA(b=None, d=1)%p50(b=None, d=1)**cyt)
Observable('o_A50_nuc', RelA(b=None, d=1)%p50(b=None, d=1)**nuc)
Observable('o_AA_cyt', RelA(b=None, d=1)%RelA(b=None, d=1)**cyt)
Observable('o_AA_nuc', RelA(b=None, d=1)%RelA(b=None, d=1)**nuc)
Observable('o_p5050_cyt', p50(b=None, d=1)%p50(b=None, d=1)**cyt)
Observable('o_p5050_nuc', p50(b=None, d=1)%p50(b=None, d=1)**nuc)
Observable('o_IkBa_cyt', IkBa(b=None, s='u')**cyt)
Observable('o_IkBa_nuc', IkBa(b=None, s='u')**nuc)
Observable('o_IkBa_p', IkBa(b=None, s='p'))
Observable('o_IkBaNFkB_cyt', (IkBa(b=1, s='u')%RelA(b=1))**cyt)
Observable('o_IkBaNFkB_nuc', (IkBa(b=1, s='u')%RelA(b=1))**nuc)
Observable('o_tIkBa', tIkBa()**cyt)
Observable('o_tA20', tA20()**cyt)
Observable('o_A20', A20()**cyt)
Observable('o_tTarget', tTarget()**cyt)
Observable('o_tComp', tp50()**nuc)



#initial conditions parameters
Parameter('IkBa_0', 374353) #chose middle range value 0.3 uM from Suzan's model and converted to molecules (0.3e-6*avagodro*cytVol)
Parameter('RelA_0', 374353)
Parameter('p50_0', 374353)
Parameter('neut_IKK_0', 100000) # .08e-6*AvN*cytv
 

#initial conditions - set RelA, p50, and IkBa to same amount
Initial(IKK(b=None, s='n')**cyt, neut_IKK_0)
Initial(IkBa(b=None, s='u')**cyt, IkBa_0)
Initial(RelA(b=None, d=None)**cyt, RelA_0)
Initial(p50(b=None, d=None)**cyt, p50_0)


#IKK and IkBa phosphorylation parameters
Parameter('kact_IKK', 0.000) # ka, set to 0.001 when TNF present, 1/s
Parameter('kinact_IKK', 0.003) #ki, 1/s
#Parameter('kneut_IKK', 0.0006) #kp, 1/s
Parameter('kphos_IkBa', 1.2288e-19) # (kc1a*1e6)/AvN, 1/M*s*AvN IKK + IkBa
Parameter('kphos_IkBaRel', 6.144e-19 ) # (kc2a*1e6)/AvN, 1/M*s*AvN IKK + IkBa:NFkB
Parameter('kdeg_pIkBa', 0.1) #kt1a and kt2a, 1/s

#dimerization parameters
#values are from Tsui Hoffmann 2011 (all Kon divided by AvN)

Parameter('kbind_RelA_p50', 5.2473e-20)
Parameter('kdis_RelA_p50', 3.16e-4)
Parameter('kbind_RelA_RelA', 1.6605e-20)
Parameter('kdis_RelA_RelA', 8e-3)
Parameter('kbind_p50_p50', 4.9816e-20)
Parameter('kdis_p50_p50', 9e-4)

#IkBa NFkB binding parameters -
#Kd from D2FC for IkBa-NFkB is 100 nM
#will keep koff the same and change kon of RelARelA and p50p50 to match ratio
#reported in Phelps et al 2000
#kd for RelARelA is 2.5 fold lower = 100 nM * 2.5  = 250 nM; kon = 2e5 (ms)**-1
#kd for p50p50 is 20 fold lower = 100 nM * 20  = 2 uM; kon = 2.5e4 (ms)**-1
#divide kon by AvN for correct units
Parameter('kbind_IkBa_A50', 8.3027e-19) # ka1a/AvN, 1/M*s*AvN IkBa + NFkB binding
Parameter('kdis_IkBa_A50', 0.05) #kd1a, 1/s IkBa:NFkB dissociation
Parameter('kbind_IkBa_AA', 3.3211e-19) # ka1a/AvN, 1/M*s*AvN IkBa + NFkB binding
Parameter('kdis_IkBa_AA', 0.05) #kd1a, 1/s IkBa:NFkB dissociation
Parameter('kbind_IkBa_p5050', 4.1513e-20) # ka1a/AvN, 1/M*s*AvN IkBa + NFkB binding
Parameter('kdis_IkBa_p5050', 0.05) #kd1a, 1/s IkBa:NFkB dissociation

#import and export
#for now, assume import and export the same for all dimers
Parameter('kimp_A50', 0.0026) #ki1, 1/s
Parameter('kexp_A50', cn*0.000052) #ke1, 1/s
Parameter('kimp_AA', 0.0026) #ki1, 1/s
Parameter('kexp_AA', cn*0.000052) #ke1, 1/s
Parameter('kimp_p5050', 0.0026) #ki1, 1/s
Parameter('kexp_p5050', cn*0.000052) #ke1, 1/s
Parameter('kexp_IkBaRel', cn*0.01) #ke2a, 1/s
Parameter('kimp_IkBa', 0.00067) #ki3a, 1/s
Parameter('kexp_IkBa', cn*0.000335) #ke3a, 1/s

#transcription and synthesis
Parameter('ksynth_IkBa', 0.5) #c2a, 1/s
Parameter('kdeg_tIkBa', 0.0003) #c3a, 1/s
Parameter('kdeg_IkBa', 0.0005) #c4a, 1/s
Parameter('kdeg_bIkBa', 0.000022) #c5a 1/s

Parameter('ksynth_A20', 0.5) #c2 1/s
Parameter('kdeg_tA20', 0.0004) #c3 1/s
Parameter('kdeg_A20', 0.0045) #c4 1/s
Parameter('kdeg_tTarget', 0.0004) #c3 1/s
Parameter('ksynth_p50', 0.5) #c2a 1/s
Parameter('kdeg_tp50', 0.000042857) #c6a 1/s
Parameter('kdeg_p50', 0.0005) #c4a 1/s

# DNA binding parameters
Parameter('dbind_A50', 24579)
Parameter('dbind_AA', 24579)
Parameter('dbind_p5050', 24579)
Parameter('dbind_Targ_p5050', 12289.5)

#replace expressions
Parameter('ktran_IkBa', 0)
Parameter('ktran_A20', 0)
Parameter('ktran_Target', 0)
Parameter('ktran_p50', 0)
Parameter('kneut_IKK', 0)

#dimerization
Rule('bind_RelA_RelA', RelA(b=None, d=None) + RelA(b=None, d=None) |
        RelA(b=None, d=1)%RelA(b=None, d=1), kbind_RelA_RelA, kdis_RelA_RelA)
Rule('bind_RelA_p50', RelA(b=None, d=None) + p50(b=None, d=None) |
        RelA(b=None, d=1)%p50(b=None, d=1), kbind_RelA_p50, kdis_RelA_p50)
Rule('bind_p50_p50', p50(b=None, d=None) + p50(b=None, d=None) |
        p50(b=None, d=1)%p50(b=None, d=1), kbind_p50_p50, kdis_p50_p50)

#IKK activation
Rule('act_IKK', IKK(b=None, s='n') >> IKK(b=None, s='a'), kact_IKK) # IKK activation
Rule('inact_IKK', IKK(b=None, s='a') >> IKK(b=None, s='i'), kinact_IKK) # IKK inactivation
Rule('neut_IKK', IKK(b=None, s='i') >> IKK(b=None, s='n'), kneut_IKK) #IKK return to neutral

#IkBa
Rule('bind_IkBa_A50', IkBa(b=None, s='u') + RelA(b=None, d=1)%p50(b=None, d=1) | 
        IkBa(b=2, s='u')%RelA(b=2, d=1)%p50(b=None, d=1), kbind_IkBa_A50, kdis_IkBa_A50) #IkBa bind RelA:p50
Rule('bind_IkBa_AA', IkBa(b=None, s='u') + RelA(b=None, d=1)%RelA(b=None, d=1) | 
        IkBa(b=2, s='u')%RelA(b=2, d=1)%RelA(b=None, d=1), kbind_IkBa_AA, kdis_IkBa_AA) #IkBa bind RelA:RelA
Rule('bind_IkBa_p5050', IkBa(b=None, s='u') + p50(b=None, d=1)%p50(b=None, d=1) | 
        IkBa(b=2, s='u')%p50(b=2, d=1)%p50(b=None, d=1), kbind_IkBa_p5050, kdis_IkBa_p5050) #IkBa bind p50:p50
Rule('phos_IkBa', IKK(b=None, s='a') + IkBa(b=None, s='u') >> 
     IKK(b=None, s='a') + IkBa(b=None, s='p'), kphos_IkBa) #phosphorylation IkBa
Rule('phos_IkBaNF', IKK(b=None, s='a') + IkBa(b=ANY, s='u') >> 
     IKK(b=None, s='a') + IkBa(b=ANY, s='p'), kphos_IkBaRel) #phosphorylation IkBa:NFkB
Rule('deg_pIkBa', IkBa(b=None, s='p') >> None, kdeg_pIkBa) #pIkBa degradation
Rule('deg_pIkBaA50', IkBa(b=2, s='p')%RelA(b=2, d=1)%p50(b=None, d=1) >>
     RelA(b=None, d=1)%p50(b=None, d=1), kdeg_pIkBa) #bound pIkBa degradation
Rule('deg_pIkBaAA', IkBa(b=2, s='p')%RelA(b=2, d=1)%RelA(b=None, d=1) >>
     RelA(b=None, d=1)%RelA(b=None, d=1), kdeg_pIkBa) #bound pIkBa degradation
Rule('deg_pIkBap5050', IkBa(b=2, s='p')%p50(b=2, d=1)%p50(b=None, d=1) >>
     p50(b=None, d=1)%p50(b=None, d=1), kdeg_pIkBa) #bound pIkBa degradation

#transport rules
Rule('import_A50', (RelA(b=None, d=1)%p50(b=None, d=1))**cyt | 
        (RelA(b=None, d=1)%p50(b=None, d=1))**nuc, kimp_A50, kexp_A50) #RelA:p50 nuclear import/export
Rule('import_AA', (RelA(b=None, d=1)%RelA(b=None, d=1))**cyt | 
        (RelA(b=None, d=1)%RelA(b=None, d=1))**nuc, kimp_AA, kexp_AA) #RelA:RelA nuclear import/export
Rule('import_p5050', (p50(b=None, d=1)%p50(b=None, d=1))**cyt | 
        (p50(b=None, d=1)%p50(b=None, d=1))**nuc, kimp_p5050, kexp_p5050) #p50:p50 nuclear import/export
Rule('export_IkBaA50', (IkBa(b=2, s='u')%RelA(b=2, d=1)%p50(b=None, d=1))**nuc >> 
     (IkBa(b=2, s='u')%RelA(b=2, d=1)%p50(b=None, d=1))**cyt, kexp_IkBaRel) #IkBa:RelA:p50 export
Rule('export_IkBaAA', (IkBa(b=2, s='u')%RelA(b=2, d=1)%RelA(b=None, d=1))**nuc >> 
     (IkBa(b=2, s='u')%RelA(b=2, d=1)%RelA(b=None, d=1))**cyt, kexp_IkBaRel) #IkBa:RelA:p50 export
Rule('export_IkBa5050', (IkBa(b=2, s='u')%p50(b=2, d=1)%p50(b=None, d=1))**nuc >> 
     (IkBa(b=2, s='u')%p50(b=2, d=1)%p50(b=None, d=1))**cyt, kexp_IkBaRel) #IkBa:RelA:p50 export
Rule('import_IkBa', IkBa(b=None, s='u')**cyt | 
        IkBa(b=None, s='u')**nuc, kimp_IkBa, kexp_IkBa) #IkBa nuclear import/export

#IkBa synthesis and degradation
Rule('tran_IkBa', None >> tIkBa()**cyt, ktran_IkBa) #IkbA mRNA transcription
Rule('synth_IkBa', tIkBa()**cyt >> tIkBa()**cyt + IkBa(b=None, s='u')**cyt, ksynth_IkBa) #IkBa protein synthesis
Rule('deg_tIkBa', tIkBa() >> None, kdeg_tIkBa) #tIkBa degradation
Rule('deg_IkBa', IkBa(b=None, s='u') >> None, kdeg_IkBa) # unbound IkBa degradation
Rule('deg_bIkBa', (IkBa(b=1, s='u')%RelA(b=1))**cyt >> RelA(b=None), kdeg_bIkBa) # bound IkBa degradation

#A20 synthesis and degradation
Rule('tran_A20', None >> tA20()**cyt, ktran_A20) #A20 mRNA transcription
Rule('synth_A20', tA20()**cyt >> tA20()**cyt + A20()**cyt, ksynth_A20) #A20 protein synthesis
Rule('deg_tA20', tA20() >> None, kdeg_tA20) #A20 transcript degradation
Rule('deg_A20', A20() >> None, kdeg_A20) #A20 protein degradation

#prototypical inducible target
Rule('tran_Target', None >> tTarget()**cyt, ktran_Target) #inducible target transcription
Rule('deg_tTarget', tTarget() >> None, kdeg_tTarget) #tTarget degradation

#Competitor protein synthesis and degradation
Rule('tran_p50', None >> tp50()**cyt, ktran_p50) #Competitor mRNA transcription
Rule('synth_p50', tp50()**cyt >> tp50()**cyt + p50(b=None, d=None)**cyt, ksynth_p50) #Comp protein synthesis
Rule('deg_tp50', tp50() >> None, kdeg_tp50) #Comp transcript degradation
Rule('deg_p50', p50(b=None, d=None) >> None, kdeg_p50) #Comp protein degradation
