# -*- coding: utf-8 -*-
"""
Created on Wed May 23 08:11:19 2018

@author: Amy Thurber
"""

from pysb import *
import numpy as np

Model()

#Create multiple species
def create_monomers(list):
	for i in list:
		x = Monomer(i,['b'])

# E + S <-> E:S -> E + Sa
def act2(enz,sub,suba,kf,kr,kc):
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
	Rule(r2_name, ES >> E + Sa, z)

# L + R <-> L:R -> Ra
def act3(enz,sub,suba,kf,kr,kc):
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
def synth_enz(enz,prod,k):
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
def migrate2(s1,s2,kf,kr):
	r_name = 'mig_%s_%s' % (s1.name, s2.name)
	x = Parameter(r_name + '_kf',kf)
	y = Parameter(r_name + '_kr',kr)
	Rule(r_name, s1(b=None) | s2(b=None), x, y)

# S1 -> S2
def migrate(s1,s2,k):
	r_name = 'mig_%s_%s' % (s1.name, s2.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, s1(b=None) >> s2(b=None), x)

Monomer('p65',['b','s','l'],{'s':['u','p'],'l':['c','n']})
Monomer('p105_p50',['s'],{'s':['u','p']})
Monomer('p50',['b','b_DNA','b_IkB','l'],{'l':['c','n']})
Monomer('p105',['b'])
Monomer('IkB',['b','s','l'],{'s':['u','p'],'l':['c','n']})
create_monomers(['DISC','TNF','TNFR','TNFRa','MKP','TAK','TAKp','TB','PAMP','TLR','TLRa','IKK','IKKp','p38','p38p','HSP27','HSP27p','HSP27n'])
create_monomers(['LXA4'])
Monomer('A20',['b'])
create_monomers(['Cox','Pol','Pol_Cox','tCox','pCox','Pol_IkB','tIkB','Pol_A20','tA20','Pol_TNF','tTNF'])
Monomer('dIkB',['b_p50','b_Pol'])
Monomer('dA20',['b_p50','b_Pol'])
Monomer('dTNF',['b_p50','b_Pol'])
Monomer('dCox',['b_p50','b_Pol'])
Monomer('dSyt',['b_p50','b_Pol'])
create_monomers(['tSyt','Syt','Pol_Syt'])
create_monomers(['PLA2','PLA2p','AA','PGE2','LO','EP4','EP4a','PIP3','AKT','AKTp','LOp','EP2','EP2a','cAMP','PKA','PKAa'])
create_monomers(['TNFR2','TNFR2a'])

Observable('o_p65',p65(b=None,s='u',l='c'))
Observable('o_p50',p50(b=None,b_DNA=None,b_IkB=None,l='c'))
Observable('o_p105',p105(b=None))
Observable('o_p65_p105',p65(b=1,s='u',l='c')%p105(b=1))
Observable('o_p65_p50',p65(b=1,s='u',l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None))
Observable('o_p65_p50n',p65(b=1,s='u',l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None))
Observable('o_p50_p50',p50(b=1,b_DNA=None,b_IkB=None,l='c')%p50(b=1,b_DNA=None,b_IkB=None,l='c'))
Observable('o_p50_p50n',p50(b=1,b_DNA=None,b_IkB=None,l='n')%p50(b=1,b_DNA=None,b_IkB=None,l='n'))
Observable('o_p65n',p65(b=None,s='u',l='n'))
Observable('o_p50n',p50(b=None,b_DNA=None,b_IkB=None,l='n'))
Observable('o_p105_p50',p105_p50(s='u'))
Observable('o_IkB',IkB(b=None,s='u',l='c'))
Observable('o_IkBn',IkB(b=None,s='u',l='n'))
Observable('o_IkB_NFkB',IkB(b=2,s='u',l='c') % p65(b=1,s='u',l='c') % p50(b=1,b_IkB=2,b_DNA=None,l='c'))
Observable('o_IkB_NFkBn',IkB(b=2,s='u',l='n') % p65(b=1,s='u',l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n'))
Observable('o_IKK',IKK(b=None))
Observable('o_TAK',TAKp(b=None))
Observable('o_TLR',TLRa(b=None))
Observable('o_A20',A20(b=None))
Observable('o_A20_IKK',A20(b=1)%IKK(b=1))
Observable('o_IKKp',IKKp(b=WILD))
Observable('o_TNF',TNF(b=WILD))
Observable('o_NFkB',p65(b=1,s='p',l='n') % p50(b=1,b_IkB=None,b_DNA=None,l='n'))
Observable('o_HSP27',HSP27n(b=WILD))
Observable('o_p38',p38p(b=WILD))
Observable('o_TNFR',TNFRa(b=WILD))
Observable('o_LXA4',LXA4(b=None))
Observable('o_AKT',AKTp(b=WILD))
Observable('o_Cox',Cox(b=WILD))
Observable('o_pCox',pCox(b=WILD))
Observable('o_tCox',tCox(b=WILD))
Observable('o_EP2',EP2a(b=WILD))
Observable('o_EP4',EP4a(b=WILD))
Observable('o_PGE2',PGE2(b=WILD))
Observable('o_PKA',PKAa(b=WILD))
Observable('o_Syt7',Syt(b=WILD))
Observable('o_DISC',DISC(b=WILD))
Observable('o_AA',AA(b=WILD))
Observable('o_tA20',tA20(b=None))
Observable('o_tIkB',tIkB(b=None))

#NF-kB module

Parameter('synth_p65_k',3*30)   #p65 synthesis rate
Parameter('kdeg_generic',3*3e-4)    #degradation rate p65,p50,p105,all dimers,IkB,IkBp65p50
Parameter('synth_p50_k',3*15)
Parameter('p105_process_k',3*3e-4) 
Parameter('kf_generic',3*1e-6) #NFkB dimer binding
Parameter('kr_generic',3*1e-3) #NFkB dimer dissociation
Parameter('mig_n_generic',3*0.0026) #nuc import
Parameter('mig_c_generic',3*0.00052) #nuc export
Parameter('mig_n_IkB',3*0.00067) #IkB nuc transport
Parameter('mig_c_IkB',3*0.000335) #IkB nuc export
Parameter('bind_IkB_kf',3*4e-7) # IkB + p65:p50 binding
Parameter('bind_IkB_kr',3*5e-4) # IkB + p65:p50 dissociation
Parameter('mig_IkB_NFkB_k',3*0.01) #IkB:p65:p50 nuc export

# initial value parameters, let model reach equilibrium w/o IKK act
Parameter('stable_p65',0)
Parameter('stable_p50',0)
Parameter('stable_p105',0)
Parameter('stable_p65_p105',0)
Parameter('stable_p65_p50',0)
Parameter('stable_p65_p50n',0)
Parameter('stable_p50_p50',0)
Parameter('stable_p50_p50n',0)
Parameter('stable_p65n',0)
Parameter('stable_p50n',0)
Parameter('stable_p105_p50',0)
Parameter('stable_IkB',0)
Parameter('stable_IkBn',0)
Parameter('stable_IkB_NFkB',0)
Parameter('stable_IkB_NFkBn',0)

Initial(p65(b=None,s='u',l='c'),stable_p65)
Initial(p50(b=None,b_DNA=None,b_IkB=None,l='c'),stable_p50)
Initial(p105(b=None),stable_p105)
Initial(p65(b=1,s='u',l='c')%p105(b=1),stable_p65_p105)
Initial(p65(b=1,s='u',l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None),stable_p65_p50)
Initial(p65(b=1,s='u',l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None),stable_p65_p50n)
Initial(p50(b=1,l='c',b_DNA=None,b_IkB=None)%p50(b=1,l='c',b_DNA=None,
        b_IkB=None),stable_p50_p50)
Initial(p50(b=1,l='n',b_DNA=None,b_IkB=None)%p50(b=1,l='n',b_DNA=None,
        b_IkB=None),stable_p50_p50n)
Initial(p65(b=None,s='u',l='n'),stable_p65)
Initial(p50(b=None,b_DNA=None,b_IkB=None,l='n'),stable_p50n)
Initial(p105_p50(s='u'),stable_p105_p50)
Initial(IkB(b=None,s='u',l='c'),stable_IkB)
Initial(IkB(b=None,s='u',l='n'),stable_IkBn)
Initial(IkB(b=1,s='u',l='c') % p65(b=2,s='u',l='c') % p50(b=2,b_IkB=1,
        b_DNA=None,l='c'),stable_IkB_NFkB)
Initial(IkB(b=1,s='u',l='n') % p65(b=2,s='u',l='n') % p50(b=2,b_IkB=1,
        b_DNA=None,l='n'),stable_IkB_NFkBn)



Rule('synth_p65', None >> p65(b=None,s='u',l='c'),synth_p65_k) #p65 synthesis
Rule('deg_p65',p65(b=None,l='c') >> None, kdeg_generic) #p65 degradatio
Rule('synth_p50',None >> p105_p50(s='u'),synth_p50_k) #p105_p50 synthesis WHAT IS THIS!!!
Rule('p105_process',p105_p50(s='u') >> p105(b=None) + p50(b=None,b_IkB=None,b_DNA=None,l='c'),
     p105_process_k) #62 p105 processing to p105 + p50
Rule('deg_p105',p105(b=None) >> None,kdeg_generic) #p105 degradation
Rule('deg_p50',p50(b=None,b_IkB=None,b_DNA=None,l='c') >> None,kdeg_generic) #p50 degradation
Rule('bind_p65_p105',p65(b=None,l='c') + p105(b=None) | p65(b=1,l='c')%p105(b=1),
     kf_generic, kr_generic) #p65+p105 binding
Rule('deg_p65_p105',p65(b=1)%p105(b=1) >> None, kdeg_generic) #p65:p105 degradation
Rule('bind_p65_p50',p65(b=None,l='c') + p50(b=None,b_IkB=None,l='c') | 
        p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c'),kf_generic,kr_generic) #p65+p50 binding
Rule('bind_p65_p50n',p65(b=None,l='n') + p50(b=None,b_DNA=None,b_IkB=None,l='n') | 
        p65(b=1,l='n')%p50(b=1,b_DNA=None,b_IkB=None,l='n'),kf_generic,kr_generic) #p65+p50 in nuc
Rule('deg_p65_p50',p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c') >> None, kdeg_generic) #p65:p50 deg
Rule('bind_p50_p50',p50(b=None,l='c') + p50(b=None,l='c') | 
        p50(b=1,l='c')%p50(b=1,l='c'),kf_generic,kr_generic) #p50+p50 binding
Rule('bind_p50_p50n',p50(b=None,b_DNA=None,l='n') + p50(b=None,b_DNA=None,l='n') | 
        p50(b=1,l='n',b_DNA=None,)%p50(b=1,l='n',b_DNA=None),kf_generic,kr_generic) #p50+p50 nuc
Rule('deg_p50_p50',p50(b=1,l='c')%p50(b=1,l='c') >> None, kdeg_generic) #p50:p50 deg
Rule('mig_p65_p50',p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None) | 
        p65(b=1,l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None),mig_n_generic,mig_c_generic) #p65:p50 nuc transport
Rule('mig_p50_p50',p50(b=1,l='c',b_DNA=None)%p50(b=1,l='c',b_DNA=None) | 
        p50(b=1,l='n',b_DNA=None)%p50(b=1,l='n',b_DNA=None),mig_n_generic,mig_c_generic) #p50:p50 nuc transport
Rule('mig_p65',p65(b=None,l='c') | p65(b=None,l='n'), mig_n_generic, mig_c_generic) #p65 nuc transport
Rule('mig_p50',p50(b=None,b_DNA=None,l='c') | 
        p50(b=None,b_DNA=None,l='n'), mig_n_generic, mig_c_generic) #p50 nuc transport
Rule('mig_IkB',IkB(b=None,s='u',l='c') | IkB(b=None,s='u',l='n'), mig_n_IkB, mig_c_IkB) #IkB nuc transport
Rule('bind_IkB_c',IkB(b=None,s='u',l='c') + p65(b=1,l='c') % p50(b=1,b_DNA=None,b_IkB=None,l='c') | 
        IkB(b=2,s='u',l='c') % p65(b=1,l='c') % p50(b=1,b_DNA=None,b_IkB=2,l='c'), bind_IkB_kf, bind_IkB_kr) #IkB bind p65:p50 cyt
Rule('bind_IkB_n',IkB(b=None,s='u',l='n') + p65(b=1,l='n') % p50(b=1,b_IkB=None,b_DNA=None,l='n') | 
        IkB(b=2,s='u',l='n') % p65(b=1,l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n'), bind_IkB_kf, bind_IkB_kr) #IkB bind p65:p50 nuc
Rule('mig_IkB_NFkB',IkB(b=2,l='n') % p65(b=1,l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n') >> 
     IkB(b=2,l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,b_DNA=None,l='c'), mig_IkB_NFkB_k) #IkB:p65:p50 nuc export
Rule('deg_IkB_NFkB',IkB(b=2,s='u',l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,l='c') >> None, kdeg_generic) #IkB:p65:p50 degradation

#IKK Reactions (called inflammatory cell signaling in thesis)

Parameter('phos_p65_kf',3e-5)
Parameter('phos_p65_kr',3e-1)
Parameter('phos_p65_kc',3e-1)
Parameter('IKKp_IkB_kf',0)
Parameter('IKKp_kf',3e6)
Parameter('inh_LXA4',1e3)
Parameter('kinase_kr_generic',3*1e-1)
Parameter('kinase_kc_generic',3*1e-1)
Parameter('IKKp_p105_kf',3*1e-6)
Parameter('p105p_process_k',3*1e-1)
Parameter('deg_IkBp_k',3*0.1)
Parameter('bind_A20_kf',3*1e-4)
Parameter('bind_A20_kr',3*1e-2)

# initial parameters
Parameter('stable_IKK',1e5)
Parameter('TAK_0',1e6)
Parameter('TLR_0',1e4)
Parameter('TB_0',1)
Parameter('LXA4_0',0)
Parameter('stable_A20',0)
Parameter('stable_A20_IKK',0)
Parameter('p38_0',1e6)
Parameter('HSP27_0',1e4)
Parameter('MKP_0',1e7)
Parameter('TNFR_0',1e4)

Initial(IKK(b=None),stable_IKK)
Initial(A20(b=None),stable_A20)
Initial(A20(b=1)%IKK(b=1),stable_A20_IKK)



synth_enz(TB,PAMP,1e2) #73 synth PAMP
deg(PAMP,1e-3) #74 degrade PAMP
act3(PAMP,TLR,TLRa,1e-6,1e-1,1e-1) #75 catalyze/activate TLRa
migrate(TLRa,TLR,1e-3) #76 inactivate TLRa
act2(TLRa,TAK,TAKp,1e-7,1e-1,1e-1) #77 activate TAK
synth(TNFR,10) #78 synthesize TNFR
deg(TNFR,1e-3) #79 degrade TNFR
migrate(TNFRa,DISC,1e-3) #82 TNFR becomes DISC
deg(DISC,1e-3) #83 degrade DISC
act3(TNF,TNFR,TNFRa,1e-7,1e-1,1e-1) #80 activate TNFR
act2(TNFRa,TAK,TAKp,0,1e-1,1e-1) #81 phos TAK
act2(MKP,TAKp,TAK,1e-6,1e-1,3e-3) #89 dephos TAK
act2(TAKp,IKK,IKKp,3*1e-6,3*1e-1,3*1e-1) #84 phos IKK
#act2(MKP,IKKp,IKK,3*1e-5,3*1e-1,3*1e-1)
act2(TAKp,p38,p38p,1e-7,3e-1,3e-1) #87 phos p38
act2(MKP,p38p,p38,3e-6,3e-1,9e-3) #90 dephos p38
#migrate(p38p,p38,1e-4)
act2(p38p,HSP27,HSP27p,3e-6,3e-1,3e-1) 
act2(MKP,HSP27p,HSP27,3e-6,3e-1,9e-3)
migrate2(HSP27p,HSP27n,0.1,0.001)

Rule('bind_p65n_HSP27n',HSP27n(b=None) + p65(l='n',s='u') | 
        HSP27n(b=1) % p65(l='n',s=('u',1)), phos_p65_kf, phos_p65_kr)
Rule('cat_p65n_HSP27n',HSP27n(b=1) % p65(s=('u',1)) >> HSP27n(b=None) + p65(s='p'),phos_p65_kc)

Expression('IKKp_IkB_kf_LXA4', 1/(IKKp_kf+inh_LXA4*o_LXA4))
Rule('bind_IKKp_IkB_LXA4', IKKp(b=None) + IkB(l='c',s='u') | 
        IKKp(b=1)%IkB(l='c',s=('u',1)), IKKp_IkB_kf_LXA4, kinase_kr_generic) #replaces rule 92, note IkB(b) is not defined)
#Rule('bind_IKKp_IkB', IKKp(b=None) + IkB(l='c',s='u') | IKKp(b=1) % 
#     IkB(l='c',s=('u',1)), IKKp_IkB_kf, kinase_kr_generic) #92 IKKp-Ikb binding (kf set to 0)
Rule('cat_IKKp_IkB', IKKp(b=1) % IkB(l='c',s=('u',1)) >> IKKp(b=None) + IkB(l='c',s='p'), kinase_kc_generic) #92
Rule('bind_IKKp_p105', IKKp(b=None) + p105_p50(s='u') | 
        IKKp(b=1) % p105_p50(s=('u',1)), IKKp_p105_kf, kinase_kr_generic) #93 IKKp-p105 binding
Rule('cat_IKKp_p105', IKKp(b=1) % p105_p50(s=('u',1)) >> IKKp(b=None) + p105_p50(s='p'), kinase_kc_generic) #93
Rule('p105p_process',p105_p50(s='p') >> p105(b=None) + p50(b=None,b_IkB=None,b_DNA=None,l='c'), p105p_process_k) #NOT IN THESIS
Rule('deg_IkBp',IkB(b=None,s='p',l='c') >> None, deg_IkBp_k) #94 degrade IkBp
Rule('deg_IkBp_bound',IkB(b=2,s='p',l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,l='c') >> 
     p65(b=1,l='c') % p50(b=1,b_IkB=None,l='c'), deg_IkBp_k) #96 p65:p50 release from IkBp
Rule('bind_A20',IKK(b=None) + A20(b=None) | IKK(b=1) % A20(b=1), bind_A20_kf, bind_A20_kr) #85 A20 IKK binding
synth(IKK,3*30) #86 synth IKK
deg(IKKp,3*3e-4) #83 degrade IKKp
deg(IKK,3*3e-4) #83 degrade IKK
Rule('deg_IKK_A20',IKK(b=1) % A20(b=1) >> None, kdeg_generic) #83 degrade IKK:A20

#Gene Expression

Parameter('Pol_kr',3*.1)
Parameter('Pol_kc',3*.1)
Parameter('transcribe_k',3*1e-1)
Parameter('translate_k',3*0.5)
Parameter('p50_DNA_kf_generic',3*1.23e-6)
Parameter('p50_DNA_kr_generic',3*1e-1)
Parameter('Pol_dIkB_kf',3*5.25e-7)
Parameter('Pol_dA20_kf',3*1e-5)
Parameter('deg_tA20_k',3*1e-3)
Parameter('deg_A20_k',3*3e-4)
Parameter('p50_TNF_kf',8e-7)
Parameter('Pol_dTNF_kf',3*1.6e-8)
Parameter('Pol_dTNF_p_kf',3*1.6e-6)
Parameter('deg_TNF_k',1e-3)
Parameter('deg_tTNF_k',0.001)
Parameter('p50_Cox_kf',8e-7)
Parameter('Pol_dCox_kf',3*1.6e-8)
Parameter('Pol_dCox_p_kf',3*1.6e-6)
Parameter('deg_Cox_k',1e-2)
Parameter('deg_tCox_k',1e-2)
Parameter('deg_pCox_k',1e-4)
Parameter('p50_Syt_kf',8e-7)
Parameter('Pol_dSyt_kf',3*1.6e-8)
Parameter('Pol_dSyt_p_kf',3*1.6e-6)
Parameter('deg_Syt_k',1e-3)
Parameter('deg_tSyt_k',5e-3)

#initial parameters
Parameter('stable_tA20',0)
Parameter('stable_tIkB',0)
Parameter('init_Pol',2e5)
Parameter('init_DNA',2)


Initial(tA20(b=None),stable_tA20)
Initial(tIkB(b=None),stable_tIkB)
Initial(Pol(b=None),init_Pol)
Initial(dIkB(b_p50=None,b_Pol=None),init_DNA)
Initial(dA20(b_p50=None,b_Pol=None),init_DNA)
Initial(dTNF(b_p50=None,b_Pol=None),init_DNA)
Initial(dCox(b_p50=None,b_Pol=None),init_DNA)
Initial(dSyt(b_p50=None,b_Pol=None),init_DNA)





#IkB transcription
Rule('bind_p50_dIkB', dIkB(b_p50=None) + p65(b=1,l='n') % p50(l='n',b=1,b_IkB=None,b_DNA=None) |
        dIkB(b_p50=2) % p65(b=1,l='n') % p50(l='n',b=1,b_IkB=None,b_DNA=2), p50_DNA_kf_generic, p50_DNA_kr_generic)
Rule('bind_Pol_dIkB',Pol(b=None) + dIkB(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dIkB(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_dIkB_kf,Pol_kr)
Rule('cat_Pol_dIkB',Pol(b=3) % dIkB(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> 
     Pol_IkB(b=None) + dIkB(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_IkB',Pol_IkB(b=None) >> Pol(b=None) + tIkB(b=None), transcribe_k)
Rule('translate_IkB',tIkB(b=None) >> tIkB(b=None) + IkB(b=None,s='u',l='c'),translate_k)
Rule('deg_tIkB',tIkB(b=None) >> None, kdeg_generic)
Rule('deg_IkB',IkB(b=None,s='u',l='c') >> None, kdeg_generic)

#A20 transcription
Rule('bind_p50_dA20', dA20(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | 
        dA20(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_DNA_kf_generic, p50_DNA_kr_generic)
Rule('bind_Pol_dA20',Pol(b=None) + dA20(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dA20(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_dA20_kf,Pol_kr)
Rule('cat_Pol_dA20',Pol(b=3) % dA20(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> 
     Pol_A20(b=None) + dA20(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_A20',Pol_A20(b=None) >> Pol(b=None) + tA20(b=None), transcribe_k)
Rule('translate_A20',tA20(b=None) >> tA20(b=None) + A20(b=None),translate_k)
Rule('deg_tA20',tA20(b=None) >> None, deg_tA20_k)
Rule('deg_A20',A20(b=None) >> None, deg_A20_k)

#TNF transcription
Rule('bind_p50_dTNF', dTNF(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | 
        dTNF(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_TNF_kf, p50_DNA_kr_generic)
Rule('bind_Pol_dTNF',Pol(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dTNF_kf,Pol_kr)
Rule('bind_Pol_dTNF_p',Pol(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dTNF_p_kf,Pol_kr)
Rule('cat_Pol_dTNF',Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> 
     Pol_TNF(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_TNF',Pol_TNF(b=None) >> Pol(b=None) + tTNF(b=None), transcribe_k)
Rule('translate_TNF',tTNF(b=None) >> tTNF(b=None) + TNF(b=None),translate_k)
Rule('deg_tTNF',tTNF(b=None) >> None, deg_tTNF_k)
Rule('deg_TNF',TNF(b=None) >> None, deg_TNF_k)

#Cox transcription
Rule('bind_p50_dCox', dCox(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | 
        dCox(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_Cox_kf, p50_DNA_kr_generic)
Rule('bind_Pol_dCox',Pol(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dCox_kf,Pol_kr)
Rule('bind_Pol_dCox_p',Pol(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dCox_p_kf,Pol_kr)
Rule('cat_Pol_dCox',Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> 
     Pol_Cox(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_Cox',Pol_Cox(b=None) >> Pol(b=None) + tCox(b=None), transcribe_k)
Rule('translate_Cox',tCox(b=None) >> tCox(b=None) + Cox(b=None),translate_k)
Rule('translate_pCox',pCox(b=None) >> pCox(b=None) + Cox(b=None),translate_k)
act2(HSP27n,tCox,pCox,1e-5,1e-1,1e-1)
Rule('deg_tCox',tCox(b=None) >> None, deg_tCox_k)
Rule('deg_pCox',pCox(b=None) >> None, deg_pCox_k)
Rule('deg_Cox',Cox(b=None) >> None, deg_Cox_k)

#Syt trasncription
Rule('bind_p50_dSyt', dSyt(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | 
        dSyt(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_Syt_kf, p50_DNA_kr_generic)
Rule('bind_Pol_dSyt',Pol(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dSyt_kf,Pol_kr)
Rule('bind_Pol_dSyt_p',Pol(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | 
        Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dSyt_p_kf,Pol_kr)
Rule('cat_Pol_dSyt',Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> 
     Pol_Syt(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_Syt',Pol_Syt(b=None) >> Pol(b=None) + tSyt(b=None), transcribe_k)
Rule('translate_Syt',tSyt(b=None) >> tSyt(b=None) + Syt(b=None),translate_k)
Rule('deg_tSyt',tSyt(b=None) >> None, deg_tSyt_k)
Rule('deg_Syt',Syt(b=None) >> None, deg_Syt_k)

#Cox Reactions

#initial parameters
Parameter('PLA2_0',1e5)
Parameter('LO_0',1e5)
Parameter('EP4_0',1e4)
Parameter('AKT_0',1e5)
Parameter('EP2_0',1e4)
Parameter('PKA_0',1e5)

act2(p38p,PLA2,PLA2p,1e-5,1e-1,1e-1)
migrate(PLA2p,PLA2,1e-3)
synth_enz(PLA2p,AA,1e-2)
deg(AA,1e-2)
act2(Cox,AA,PGE2,1e-6,1e-1,1e-1)
deg(PGE2,1e-2)
act2(LO,AA,LXA4,1e-8,1e-1,1e-1)
deg(LXA4,1e-2)
act3(PGE2,EP4,EP4a,1e-6,1e-1,1e-1)
migrate(EP4a,EP4,1e-3)
synth_enz(EP4a,PIP3,1e-2)
deg(PIP3,1e-1)
act2(PIP3,AKT,AKTp,1e-6,1e-1,1e-1)
migrate(AKTp,AKT,1e-3)
act2(AKTp,IKK,IKKp,0,1e-1,1e-1)
act3(PGE2,EP2,EP2a,1e-7,1e-1,1e-1)
migrate(EP2a,EP2,1e-3)
synth_enz(EP2a,cAMP,1e-2)
deg(cAMP,1e-1)
act3(cAMP,PKA,PKAa,1e-6,1e-1,1e-1)
migrate(PKAa,PKA,1e-3)
act2(PKAa,LO,LOp,1e-6,1e-1,1e-1)
migrate(LOp,LO,1e-3)


synth(TNFR2,0)
deg(TNFR2,1e-3)
deg(TNFR2a,1e-3)
act3(TNF,TNFR2,TNFR2a,1e-7,1e-1,1e-1)
synth(TNF,0)

for m in model.monomers:
	ic_param = model.parameters.get('%s_0' % m.name)
	if ic_param is not None:
		sites = {}
		for s in m.sites:
			if s in m.site_states:
				sites[s] = m.site_states[s][0]
			else:
				sites[s] = None
		Initial(m(sites), ic_param)