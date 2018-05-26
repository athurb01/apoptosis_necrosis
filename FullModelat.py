from pysb import *
import numpy as np

Model()

###Initial Conditions###
Parameter('p38_0',1e6)
Parameter('HSP27_0',1e4)
Parameter('MKP_0',1e7)

### Helper Functions ###

#Create multiple species
def create_monomers(list):
	for i in list:
		x = Monomer(i,['b'])


#Degrade a group of monomers with a common degradation rate
def group_deg(list,k_deg):
	for m in list:
		deg(m,k_deg)

#Degrade a group of bound monomer pairs with a common degradation rate. "List" elements are in the format (monomer[0],monomer[1])
def group_degB(list,k_deg):
	for m in list:
		deg_bound(m[0],m[1],k_deg)

#Create a group of synthesis reactions for monomers with initial conditions and degradation reactions (which must already be declared):
def group_synth(list):
	for m in list:
		k_deg = model.parameters.get('deg_' + m.name + '_k').value
		IC = model.parameters.get(m.name + '_0').value
		synth(m,k_deg*IC)


# E + S -> E + Sa
def act1(enz,sub,suba,k):
	r_name = '%s_%s_cat' % (enz.name, sub.name)
	x = Parameter(r_name + '_k', k)
	Rule(r_name, enz(b=None) + sub(b=None) >> enz(b=None) + suba(b=None), x)

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

# E + S <-> E:S; E:S + ATP -> E + Sa + ADP
#(Like Act2, but ATP-dependent)
def act4(enz,sub,suba,kf,kr,kc):
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
	Rule(r2_name, ES + ATP(b=None) >> E + Sa + ADP(b=None), z)


# S -> 0
def deg(sub,k):
	r_name = 'deg_%s' % (sub.name)
	x = Parameter(r_name + '_k', k)
	Rule(r_name, sub(b=None) >> None, x)

# S1:S2 -> 0
def deg_bound(S1,S2,k):
	r_name = 'deg_%s_%s' % (S1.name,S2.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, S1(b=1)%S2(b=1) >> None, x)

# S1:S2 -> S2
def deg_bound2(S1,S2,k):
	r_name = 'deg_%s_%s' % (S1.name,S2.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, S1(b=1)%S2(b=1) >> S2(b=None), x)

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
	Rule(r_name, s1(b=None) + s2(b=None) | s3(b=None), x, y)


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


### Reaction List ###



#NF-kB module
Monomer('p65',['b','s','l'],{'s':['u','p'],'l':['c','n']})
Parameter('stable_p65',0)
Initial(p65(b=None,s='u',l='c'),stable_p65)
Observable('o_p65',p65(b=None,s='u',l='c'))
Parameter('synth_p65_k',3*30)
Rule('synth_p65', None >> p65(b=None,s='u',l='c'),synth_p65_k)
Parameter('kdeg_generic',3*3e-4)
Rule('deg_p65',p65(b=None,l='c') >> None, kdeg_generic)
Monomer('p105_p50',['s'],{'s':['u','p']})
Parameter('synth_p50_k',3*15)
Rule('synth_p50',None >> p105_p50(s='u'),synth_p50_k)
Parameter('p105_process_k',3*3e-4)
Monomer('p50',['b','b_DNA','b_IkB','l'],{'l':['c','n']})
Parameter('stable_p50',0)
Initial(p50(b=None,b_DNA=None,b_IkB=None,l='c'),stable_p50)
Observable('o_p50',p50(b=None,b_DNA=None,b_IkB=None,l='c'))
Monomer('p105',['b'])
Parameter('stable_p105',0)
Initial(p105(b=None),stable_p105)
Observable('o_p105',p105(b=None))
Rule('p105_process',p105_p50(s='u') >> p105(b=None) + p50(b=None,b_IkB=None,b_DNA=None,l='c'), p105_process_k)
Rule('deg_p105',p105(b=None) >> None,kdeg_generic)
Rule('deg_p50',p50(b=None,b_IkB=None,b_DNA=None,l='c') >> None,kdeg_generic)
Parameter('kf_generic',3*1e-6)
Parameter('kr_generic',3*1e-3)
Rule('bind_p65_p105',p65(b=None,l='c') + p105(b=None) | p65(b=1,l='c')%p105(b=1), kf_generic, kr_generic)
Parameter('stable_p65_p105',0)
Initial(p65(b=1,s='u',l='c')%p105(b=1),stable_p65_p105)
Observable('o_p65_p105',p65(b=1,s='u',l='c')%p105(b=1))
Rule('deg_p65_p105',p65(b=1)%p105(b=1) >> None, kdeg_generic)
Rule('bind_p65_p50',p65(b=None,l='c') + p50(b=None,b_IkB=None,l='c') | p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c'),kf_generic,kr_generic)
Parameter('stable_p65_p50',0)
Initial(p65(b=1,s='u',l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None),stable_p65_p50)
Observable('o_p65_p50',p65(b=1,s='u',l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None))
Rule('bind_p65_p50n',p65(b=None,l='n') + p50(b=None,b_DNA=None,b_IkB=None,l='n') | p65(b=1,l='n')%p50(b=1,b_DNA=None,b_IkB=None,l='n'),kf_generic,kr_generic)
Parameter('stable_p65_p50n',0)
Initial(p65(b=1,s='u',l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None),stable_p65_p50n)
Observable('o_p65_p50n',p65(b=1,s='u',l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None))
Rule('deg_p65_p50',p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c') >> None, kdeg_generic)
Rule('bind_p50_p50',p50(b=None,l='c') + p50(b=None,l='c') | p50(b=1,l='c')%p50(b=1,l='c'),kf_generic,kr_generic)
Parameter('stable_p50_p50',0)
Initial(p50(b=1,l='c',b_DNA=None,b_IkB=None)%p50(b=1,l='c',b_DNA=None,b_IkB=None),stable_p50_p50)
Observable('o_p50_p50',p50(b=1,b_DNA=None,b_IkB=None,l='c')%p50(b=1,b_DNA=None,b_IkB=None,l='c'))
Rule('bind_p50_p50n',p50(b=None,b_DNA=None,l='n') + p50(b=None,b_DNA=None,l='n') | p50(b=1,l='n',b_DNA=None,)%p50(b=1,l='n',b_DNA=None),kf_generic,kr_generic)
Parameter('stable_p50_p50n',0)
Initial(p50(b=1,l='n',b_DNA=None,b_IkB=None)%p50(b=1,l='n',b_DNA=None,b_IkB=None),stable_p50_p50n)
Observable('o_p50_p50n',p50(b=1,b_DNA=None,b_IkB=None,l='n')%p50(b=1,b_DNA=None,b_IkB=None,l='n'))
Rule('deg_p50_p50',p50(b=1,l='c')%p50(b=1,l='c') >> None, kdeg_generic)
Parameter('mig_n_generic',3*0.0026)
Parameter('mig_c_generic',3*0.00052)
Rule('mig_p65_p50',p65(b=1,l='c')%p50(b=1,b_IkB=None,l='c',b_DNA=None) | p65(b=1,l='n')%p50(b=1,b_IkB=None,l='n',b_DNA=None),mig_n_generic,mig_c_generic)
Rule('mig_p50_p50',p50(b=1,l='c',b_DNA=None)%p50(b=1,l='c',b_DNA=None) | p50(b=1,l='n',b_DNA=None)%p50(b=1,l='n',b_DNA=None),mig_n_generic,mig_c_generic)
Rule('mig_p65',p65(b=None,l='c') | p65(b=None,l='n'), mig_n_generic, mig_c_generic)
Parameter('stable_p65n',0)
Initial(p65(b=None,s='u',l='n'),stable_p65)
Observable('o_p65n',p65(b=None,s='u',l='n'))
Rule('mig_p50',p50(b=None,b_DNA=None,l='c') | p50(b=None,b_DNA=None,l='n'), mig_n_generic, mig_c_generic)
Parameter('stable_p50n',0)
Initial(p50(b=None,b_DNA=None,b_IkB=None,l='n'),stable_p50n)
Observable('o_p50n',p50(b=None,b_DNA=None,b_IkB=None,l='n'))
Monomer('IkB',['b','s','l'],{'s':['u','p'],'l':['c','n']})
Parameter('stable_IkB',0)
Initial(IkB(b=None,s='u',l='c'),stable_IkB)
Observable('o_IkB',IkB(b=None,s='u',l='c'))
Parameter('mig_n_IkB',3*0.00067)
Parameter('mig_c_IkB',3*0.000335)
Rule('mig_IkB',IkB(b=None,s='u',l='c') | IkB(b=None,s='u',l='n'), mig_n_IkB, mig_c_IkB)
Parameter('stable_IkBn',0)
Initial(IkB(b=None,s='u',l='n'),stable_IkBn)
Observable('o_IkBn',IkB(b=None,s='u',l='n'))
Parameter('bind_IkB_kf',3*4e-7)
Parameter('bind_IkB_kr',3*5e-4)
Rule('bind_IkB_c',IkB(b=None,s='u',l='c') + p65(b=1,l='c') % p50(b=1,b_DNA=None,b_IkB=None,l='c') | IkB(b=2,s='u',l='c') % p65(b=1,l='c') % p50(b=1,b_DNA=None,b_IkB=2,l='c'), bind_IkB_kf, bind_IkB_kr)
Parameter('stable_IkB_NFkB',0)
Initial(IkB(b=1,s='u',l='c') % p65(b=2,s='u',l='c') % p50(b=2,b_IkB=1,b_DNA=None,l='c'),stable_IkB_NFkB)
Observable('o_IkB_NFkB',IkB(b=2,s='u',l='c') % p65(b=1,s='u',l='c') % p50(b=1,b_IkB=2,b_DNA=None,l='c'))
Rule('bind_IkB_n',IkB(b=None,s='u',l='n') + p65(b=1,l='n') % p50(b=1,b_IkB=None,b_DNA=None,l='n') | IkB(b=2,s='u',l='n') % p65(b=1,l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n'), bind_IkB_kf, bind_IkB_kr)
Parameter('stable_IkB_NFkBn',0)
Initial(IkB(b=1,s='u',l='n') % p65(b=2,s='u',l='n') % p50(b=2,b_IkB=1,b_DNA=None
,l='n'),stable_IkB_NFkBn)
Observable('o_IkB_NFkBn',IkB(b=2,s='u',l='n') % p65(b=1,s='u',l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n'))
Parameter('mig_IkB_NFkB_k',3*0.01)
Rule('mig_IkB_NFkB',IkB(b=2,l='n') % p65(b=1,l='n') % p50(b=1,b_IkB=2,b_DNA=None,l='n') >> IkB(b=2,l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,b_DNA=None,l='c'), mig_IkB_NFkB_k)
Rule('deg_IkB_NFkB',IkB(b=2,s='u',l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,l='c') >> None, kdeg_generic) 



#IKK Reactions
create_monomers(['DISC','TNF','TNFR','TNFRa','MKP','TAK','TAKp','TB','PAMP','TLR','TLRa','IKK','IKKp','p38','p38p','HSP27','HSP27p','HSP27n'])
Parameter('stable_IKK',1e5)
Initial(IKK(b=None),stable_IKK)
Observable('o_IKK',IKK(b=None))
Parameter('TAK_0',1e6)
Parameter('TLR_0',1e4)
Observable('o_TAK',TAKp(b=None))
Observable('o_TLR',TLRa(b=None))
Parameter('TB_0',1)
synth_enz(TB,PAMP,1e2)
deg(PAMP,1e-3)
act3(PAMP,TLR,TLRa,1e-6,1e-1,1e-1)
migrate(TLRa,TLR,1e-3)
act2(TLRa,TAK,TAKp,1e-7,1e-1,1e-1)
synth(TNFR,10)
deg(TNFR,1e-3)
Parameter('TNFR_0',1e4)
migrate(TNFRa,DISC,1e-3)
deg(DISC,1e-3)
act3(TNF,TNFR,TNFRa,1e-7,1e-1,1e-1)
act2(TNFRa,TAK,TAKp,1e-7,1e-1,1e-1)
act2(MKP,TAKp,TAK,1e-6,1e-1,3e-3)
act2(TAKp,IKK,IKKp,3*1e-6,3*1e-1,3*1e-1)
#act2(MKP,IKKp,IKK,3*1e-5,3*1e-1,3*1e-1)
act2(TAKp,p38,p38p,1e-7,3e-1,3e-1)
act2(MKP,p38p,p38,3e-6,3e-1,9e-3)
#migrate(p38p,p38,1e-4)
act2(p38p,HSP27,HSP27p,3e-6,3e-1,3e-1)
act2(MKP,HSP27p,HSP27,3e-6,3e-1,9e-3)
migrate2(HSP27p,HSP27n,0.1,0.001)
Parameter('phos_p65_kf',3e-5)
Parameter('phos_p65_kr',3e-1)
Parameter('phos_p65_kc',3e-1)
Rule('bind_p65n_HSP27n',HSP27n(b=None) + p65(l='n',s='u') | HSP27n(b=1) % p65(l='n',s=('u',1)), phos_p65_kf, phos_p65_kr)
Rule('cat_p65n_HSP27n',HSP27n(b=1) % p65(s=('u',1)) >> HSP27n(b=None) + p65(s='p'),phos_p65_kc)
Parameter('IKKp_IkB_kf',0)
Parameter('IKKp_kf',3e6)
Parameter('inh_LXA4',1e3)
create_monomers(['LXA4'])
Parameter('LXA4_0',0)
Observable('o_LXA4',LXA4(b=None))
Parameter('kinase_kr_generic',3*1e-1)
Parameter('kinase_kc_generic',3*1e-1)
# =============================================================================
# cp = []
# cp.append((IKKp(b=1)%IkB(b=None,l='c',s=('u',1)),IkB(b=None,l='c',s='u')))
# cp.append((IKKp(b=1) % IkB(b=2, s=('u', 1), l='c') % p50(b=3, b_DNA=None, b_IkB=2, l='c') % p65(b=3, s='u', l='c'), IkB(b=1, s='u', l='c') % p50(b=2, b_DNA=None, b_IkB=1, l='c') % p65(b=2, s='u', l='c')))
# cp.append((IKKp(b=1) % IkB(b=2, s=('u', 1), l='c') % p50(b=3, b_DNA=None, b_IkB=2, l='c') % p65(b=3, s='p', l='c'), IkB(b=1, s='u', l='c') % p50(b=2, b_DNA=None, b_IkB=1, l='c') % p65(b=2, s='p', l='c')))
# 
# 
# subs=0
# for w in cp:
# 	subs=subs+1
# 	Extra(w[0],'x*y/(a+b*z)',{IKKp(b=None):'x',w[1]:'y',LXA4(b=None):'z'},{IKKp_kf:'a',inh_LXA4:'b'},'IKKp_1_' + str(subs)) 
# 	Extra(IKKp(b=None),'-x*y/(a+b*z)',{IKKp(b=None):'x',w[1]:'y',LXA4(b=None):'z'},{IKKp_kf:'a',inh_LXA4:'b'},'IKKp_2_' + str(subs))
# 	Extra(w[1],'-x*y/(a+b*z)',{IKKp(b=None):'x',w[1]:'y',LXA4(b=None):'z'},{IKKp_kf:'a',inh_LXA4:'b'},'IKKp_3_' + str(subs))
# 
# =============================================================================
# replace Robert's Extra with Expression and rule same as Frontside
Expression('IKKp_IkB_kf_LXA4', 1/(IKKp_kf+inh_LXA4*o_LXA4))
Rule('bind_IKKp_IkB_LXA4', IKKp(b=None) + IkB(l='c',s='u') | 
        IKKp(b=1)%IkB(l='c',s=('u',1)), IKKp_IkB_kf_LXA4, kinase_kr_generic) #replaces rule 92, note IkB(b) is not defined)



#Rule('bind_IKKp_IkB', IKKp(b=None) + IkB(l='c',s='u') <> IKKp(b=1) % IkB(l='c',s=('u',1)), IKKp_IkB_kf, kinase_kr_generic)
Rule('cat_IKKp_IkB', IKKp(b=1) % IkB(l='c',s=('u',1)) >> IKKp(b=None) + IkB(l='c',s='p'), kinase_kc_generic)
Parameter('IKKp_p105_kf',3*1e-6)
Rule('bind_IKKp_p105', IKKp(b=None) + p105_p50(s='u') | IKKp(b=1) % p105_p50(s=('u',1)), IKKp_p105_kf, kinase_kr_generic)
Rule('cat_IKKp_p105', IKKp(b=1) % p105_p50(s=('u',1)) >> IKKp(b=None) + p105_p50(s='p'), kinase_kc_generic)
Parameter('p105p_process_k',3*1e-1)
Rule('p105p_process',p105_p50(s='p') >> p105(b=None) + p50(b=None,b_IkB=None,b_DNA=None,l='c'), p105p_process_k)
Parameter('deg_IkBp_k',3*0.1)
Rule('deg_IkBp',IkB(b=None,s='p',l='c') >> None, deg_IkBp_k)
Rule('deg_IkBp_bound',IkB(b=2,s='p',l='c') % p65(b=1,l='c') % p50(b=1,b_IkB=2,l='c') >> p65(b=1,l='c') % p50(b=1,b_IkB=None,l='c'), deg_IkBp_k)
Monomer('A20',['b'])
Parameter('stable_A20',0)
Initial(A20(b=None),stable_A20)
Observable('o_A20',A20(b=None))
Parameter('bind_A20_kf',3*1e-4)
Parameter('bind_A20_kr',3*1e-2)
Rule('bind_A20',IKK(b=None) + A20(b=None) | IKK(b=1) % A20(b=1), bind_A20_kf, bind_A20_kr)
Parameter('stable_A20_IKK',0)
Initial(A20(b=1)%IKK(b=1),stable_A20_IKK)
Observable('o_A20_IKK',A20(b=1)%IKK(b=1))
synth(IKK,3*30)
deg(IKKp,3*3e-4)
deg(IKK,3*3e-4)
Rule('deg_IKK_A20',IKK(b=1) % A20(b=1) >> None, kdeg_generic)

#Gene Expression
create_monomers(['Cox','Pol','Pol_Cox','tCox','pCox','Pol_IkB','tIkB','Pol_A20','tA20','Pol_TNF','tTNF'])
Parameter('stable_tA20',0)
Initial(tA20(b=None),stable_tA20)
Observable('o_tA20',tA20(b=None))
Parameter('stable_tIkB',0)
Initial(tIkB(b=None),stable_tIkB)
Observable('o_tIkB',tIkB(b=None))


Parameter('init_Pol',2e5)
Initial(Pol(b=None),init_Pol)
Parameter('Pol_kr',3*.1)
Parameter('Pol_kc',3*.1)
Parameter('transcribe_k',3*1e-1)
Parameter('translate_k',3*0.5)
Parameter('init_DNA',2)

Monomer('dIkB',['b_p50','b_Pol'])
Initial(dIkB(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_DNA_kf_generic',3*1.23e-6)
Parameter('p50_DNA_kr_generic',3*1e-1)
Rule('bind_p50_dIkB', dIkB(b_p50=None) + p65(b=1,l='n') % p50(l='n',b=1,b_IkB=None,b_DNA=None) | dIkB(b_p50=2) % p65(b=1,l='n') % p50(l='n',b=1,b_IkB=None,b_DNA=2), p50_DNA_kf_generic, p50_DNA_kr_generic)
Parameter('Pol_dIkB_kf',3*5.25e-7)
Rule('bind_Pol_dIkB',Pol(b=None) + dIkB(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) | Pol(b=3) % dIkB(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_dIkB_kf,Pol_kr)
Rule('cat_Pol_dIkB',Pol(b=3) % dIkB(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_IkB(b=None) + dIkB(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_IkB',Pol_IkB(b=None) >> Pol(b=None) + tIkB(b=None), transcribe_k)
Rule('translate_IkB',tIkB(b=None) >> tIkB(b=None) + IkB(b=None,s='u',l='c'),translate_k)
Rule('deg_tIkB',tIkB(b=None) >> None, kdeg_generic)
Rule('deg_IkB',IkB(b=None,s='u',l='c') >> None, kdeg_generic)

Monomer('dA20',['b_p50','b_Pol'])
Initial(dA20(b_p50=None,b_Pol=None),init_DNA)
Rule('bind_p50_dA20', dA20(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dA20(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_DNA_kf_generic, p50_DNA_kr_generic)
Parameter('Pol_dA20_kf',3*1e-5)
Rule('bind_Pol_dA20',Pol(b=None) + dA20(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) | Pol(b=3) % dA20(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_dA20_kf,Pol_kr)
Rule('cat_Pol_dA20',Pol(b=3) % dA20(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_A20(b=None) + dA20(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_A20',Pol_A20(b=None) >> Pol(b=None) + tA20(b=None), transcribe_k)
Rule('translate_A20',tA20(b=None) >> tA20(b=None) + A20(b=None),translate_k)
Parameter('deg_tA20_k',3*1e-3)
Parameter('deg_A20_k',3*3e-4)
Rule('deg_tA20',tA20(b=None) >> None, deg_tA20_k)
Rule('deg_A20',A20(b=None) >> None, deg_A20_k)


Monomer('dTNF',['b_p50','b_Pol'])
Initial(dTNF(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_TNF_kf',8e-7)
Rule('bind_p50_dTNF', dTNF(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dTNF(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_TNF_kf, p50_DNA_kr_generic)
Parameter('Pol_dTNF_kf',3*1.6e-8)
Rule('bind_Pol_dTNF',Pol(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dTNF_kf,Pol_kr)
Parameter('Pol_dTNF_p_kf',3*1.6e-6)
Rule('bind_Pol_dTNF_p',Pol(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dTNF_p_kf,Pol_kr)
Rule('cat_Pol_dTNF',Pol(b=3) % dTNF(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_TNF(b=None) + dTNF(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_TNF',Pol_TNF(b=None) >> Pol(b=None) + tTNF(b=None), transcribe_k)
Rule('translate_TNF',tTNF(b=None) >> tTNF(b=None) + TNF(b=None),translate_k)
Parameter('deg_TNF_k',1e-3)
Parameter('deg_tTNF_k',0.001)
Rule('deg_tTNF',tTNF(b=None) >> None, deg_tTNF_k)
Rule('deg_TNF',TNF(b=None) >> None, deg_TNF_k)


Monomer('dCox',['b_p50','b_Pol'])
Initial(dCox(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_Cox_kf',8e-7)
Rule('bind_p50_dCox', dCox(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dCox(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_Cox_kf, p50_DNA_kr_generic)
Parameter('Pol_dCox_kf',3*1.6e-8)
Rule('bind_Pol_dCox',Pol(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dCox_kf,Pol_kr)
Parameter('Pol_dCox_p_kf',3*1.6e-6)
Rule('bind_Pol_dCox_p',Pol(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dCox_p_kf,Pol_kr)
Rule('cat_Pol_dCox',Pol(b=3) % dCox(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_Cox(b=None) + dCox(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_Cox',Pol_Cox(b=None) >> Pol(b=None) + tCox(b=None), transcribe_k)
Rule('translate_Cox',tCox(b=None) >> tCox(b=None) + Cox(b=None),translate_k)
Parameter('deg_Cox_k',1e-2)
Parameter('deg_tCox_k',1e-2)
Parameter('deg_pCox_k',1e-4)
Rule('translate_pCox',pCox(b=None) >> pCox(b=None) + Cox(b=None),translate_k)
act2(HSP27n,tCox,pCox,1e-5,1e-1,1e-1)
Rule('deg_tCox',tCox(b=None) >> None, deg_tCox_k)
Rule('deg_pCox',pCox(b=None) >> None, deg_pCox_k)
Rule('deg_Cox',Cox(b=None) >> None, deg_Cox_k)

Monomer('dSyt',['b_p50','b_Pol'])
create_monomers(['tSyt','Syt','Pol_Syt'])
Initial(dSyt(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_Syt_kf',8e-7)
Rule('bind_p50_dSyt', dSyt(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dSyt(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_Syt_kf, p50_DNA_kr_generic)
Parameter('Pol_dSyt_kf',3*1.6e-8)
Rule('bind_Pol_dSyt',Pol(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) | Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dSyt_kf,Pol_kr)
Parameter('Pol_dSyt_p_kf',3*1.6e-6)
Rule('bind_Pol_dSyt_p',Pol(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) | Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dSyt_p_kf,Pol_kr)
Rule('cat_Pol_dSyt',Pol(b=3) % dSyt(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_Syt(b=None) + dSyt(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_Syt',Pol_Syt(b=None) >> Pol(b=None) + tSyt(b=None), transcribe_k)
Rule('translate_Syt',tSyt(b=None) >> tSyt(b=None) + Syt(b=None),translate_k)
Parameter('deg_Syt_k',1e-3)
Parameter('deg_tSyt_k',5e-3)
Rule('deg_tSyt',tSyt(b=None) >> None, deg_tSyt_k)
Rule('deg_Syt',Syt(b=None) >> None, deg_Syt_k)

Monomer('dFLIP',['b_p50','b_Pol'])
create_monomers(['tFLIP','FLIP','Pol_FLIP'])
Initial(dFLIP(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_FLIP_kf',8e-7)
Rule('bind_p50_dFLIP', dFLIP(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dFLIP(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_FLIP_kf, p50_DNA_kr_generic)
Parameter('Pol_dFLIP_kf',0)
Rule('bind_Pol_dFLIP',Pol(b=None) + dFLIP(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) >> Pol(b=3) % dFLIP(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dFLIP_kf)
Parameter('Pol_dFLIP_p_kf',0)
Rule('bind_Pol_dFLIP_p',Pol(b=None) + dFLIP(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) >> Pol(b=3) % dFLIP(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dFLIP_p_kf)
Parameter('Pol_dFLIP_basal_kf',0)
Rule('bind_Pol_dFLIP_basal',Pol(b=None) + dFLIP(b_Pol=None,b_p50=WILD) | Pol(b=1) % dFLIP(b_Pol=1,b_p50=WILD), Pol_dFLIP_basal_kf,Pol_kr)
Rule('cat_Pol_dFLIP',Pol(b=3) % dFLIP(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_FLIP(b=None) + dFLIP(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_FLIP',Pol_FLIP(b=None) >> Pol(b=None) + tFLIP(b=None), transcribe_k)
Rule('translate_FLIP',tFLIP(b=None) >> tFLIP(b=None) + FLIP(b=None),translate_k)
synth(FLIP,5.8e-3)
Parameter('FLIP_0',2e3)
Parameter('deg_FLIP_k',2.9e-6)
Parameter('deg_tFLIP_k',1e-3)
Rule('deg_tFLIP',tFLIP(b=None) >> None, deg_tFLIP_k)
Rule('deg_FLIP',FLIP(b=None) >> None, deg_FLIP_k)

Monomer('dBcl2',['b_p50','b_Pol'])
create_monomers(['tBcl2','Bcl2','Pol_Bcl2'])
Initial(dBcl2(b_p50=None,b_Pol=None),init_DNA)
Parameter('p50_Bcl2_kf',8e-7)
Rule('bind_p50_dBcl2', dBcl2(b_p50=None) + p50(l='n',b=ANY,b_IkB=None,b_DNA=None) | dBcl2(b_p50=1) % p50(l='n',b=ANY,b_IkB=None,b_DNA=1), p50_Bcl2_kf, p50_DNA_kr_generic)
Parameter('Pol_dBcl2_kf',1e-8)
Rule('bind_Pol_dBcl2',Pol(b=None) + dBcl2(b_Pol=None,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2) >> Pol(b=3) % dBcl2(b_Pol=3,b_p50=2) % p65(b=1,s='u') % p50(b=1,b_DNA=2), Pol_dBcl2_kf)
Parameter('Pol_dBcl2_p_kf',1e-8)
Rule('bind_Pol_dBcl2_p',Pol(b=None) + dBcl2(b_Pol=None,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2) >> Pol(b=3) % dBcl2(b_Pol=3,b_p50=2) % p65(b=1,s='p') % p50(b=1,b_DNA=2), Pol_dBcl2_p_kf)
Parameter('Pol_dBcl2_basal_kf',0)
Rule('bind_Pol_dBcl2_basal',Pol(b=None) + dBcl2(b_Pol=None,b_p50=WILD) | Pol(b=1) % dBcl2(b_Pol=1,b_p50=WILD), Pol_dBcl2_basal_kf,Pol_kr)
Rule('cat_Pol_dBcl2',Pol(b=3) % dBcl2(b_Pol=3,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2) >> Pol_Bcl2(b=None) + dBcl2(b_Pol=None,b_p50=2) % p65(b=1) % p50(b=1,b_DNA=2), Pol_kc)
Rule('transcribe_Bcl2',Pol_Bcl2(b=None) >> Pol(b=None) + tBcl2(b=None), transcribe_k)
Rule('translate_Bcl2',tBcl2(b=None) >> tBcl2(b=None) + Bcl2(b=None),translate_k)
synth(Bcl2,1)
Parameter('Bcl2_0',1e4) #2e4
Parameter('deg_Bcl2_k',1e-4)
Parameter('deg_tBcl2_k',1e-3)
Rule('deg_tBcl2',tBcl2(b=None) >> None, deg_Bcl2_k)
Rule('deg_Bcl2',Bcl2(b=None) >> None, deg_Bcl2_k)



#Cox Reactions
create_monomers(['PLA2','PLA2p','AA','PGE2','LO','EP4','EP4a','PIP3','AKT','AKTp','LOp','EP2','EP2a','cAMP','PKA','PKAa'])
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
act2(LO,AA,LXA4,1e-5,1e-1,1e-1)
deg(LXA4,1e-2)
act3(PGE2,EP4,EP4a,1e-6,1e-1,1e-1)
migrate(EP4a,EP4,1e-3)
synth_enz(EP4a,PIP3,1e-2)
deg(PIP3,1e-1)
act2(PIP3,AKT,AKTp,1e-6,1e-1,1e-1)
migrate(AKTp,AKT,1e-3)
act2(AKTp,IKK,IKKp,1e-5,1e-1,1e-1)
act3(PGE2,EP2,EP2a,1e-7,1e-1,1e-1)
migrate(EP2a,EP2,1e-3)
synth_enz(EP2a,cAMP,1e-2)
deg(cAMP,1e-1)
act3(cAMP,PKA,PKAa,1e-6,1e-1,1e-1)
migrate(PKAa,PKA,1e-3)
act2(PKAa,LO,LOp,1e-6,1e-1,1e-1)
migrate(LOp,LO,1e-3)

#EARM
create_monomers(['XIAP','BAR','Mcl','C8','C8a','C3','C3a','C3i','PARP','cPARP','C6','C6a','Bid','tBid','Bax','aBax','MBax','Bax2','Bax4','MOM','MOMP','mCytoC','aCytoC','cCytoC','mSmac','aSmac','cSmac','Apaf','aApaf','C9','Apop','ATP','ADP'])
mito_vol = 0.01
Parameter('XIAP_0',1e5)
#Parameter('stable_Bcl2',0)
#Initial(Bcl2(b=None), stable_Bcl2)
#Parameter('stable_FLIP',0)
#Observable('o_Bcl2',Bcl2(b=None))
#Initial(FLIP(b=None), stable_FLIP)
#Observable('o_FLIP',FLIP(b=None))
#Parameter('stable_tFLIP',0)
#Initial(tFLIP(b=None), stable_tFLIP)
#Observable('o_tFLIP',tFLIP(b=None))
#Parameter('stable_tBcl2',0)
#Initial(tBcl2(b=None), stable_tBcl2)
#Observable('o_tBcl2',tBcl2(b=None))
Parameter('C8_0',1e4) # 2e4
Parameter('BAR_0',1e3)
Parameter('C3_0',1e4)
Parameter('C6_0',1e4)
Parameter('PARP_0',1e6)
Parameter('mCytoC_0',5e4)
Parameter('Apaf_0',1e5)
Parameter('C9_0', 1e5)
Parameter('mSmac_0',1e5)
Parameter('Bid_0',6e4) # 4e4
Parameter('Mcl_0',2e4)
Parameter('Bax_0',8e4) # 1e5
Parameter('MOM_0',5e5)
Parameter('Akt_0',2e4) 
bind(FLIP,DISC,1e-6,1e-3)
act2(DISC,C8,C8a,1e-8,1e-3,1e-5)
bind(C8a,BAR,1e-6,1e-3)
act4(C8a,C3,C3a,1e-7,1e-3,1e-8)
act4(C3a,C6,C6a,1e-7,1e-3,1e-8)
act4(C6a,C8,C8a,1e-7,1e-3,1e-8)
act4(C3a,PARP,cPARP,1e-6,1e-3,1e-8)
act2(XIAP,C3a,C3i,2e-6,1e-3,1e-1)
act4(C8a,Bid,tBid,1e-7,1e-3,1e-8)
bind(tBid,Mcl,1e-6,1e-3)
act2(tBid,Bax,aBax,1e-7,1e-3,1)
migrate2(aBax,MBax,1,1e-2)
bind(Bcl2,MBax,1e-6/mito_vol,1e-3)
bind_consume(MBax,MBax,Bax2,2e-6/mito_vol,1e-3)
bind_consume(Bax2,Bax2,Bax4,2e-6/mito_vol,1e-3)
bind(Bcl2,Bax2,1e-6/mito_vol,1e-3)
Rule('Bax2_reset',Bcl2(b=1)%Bax2(b=1) >> Bax2(b=None),deg_Bcl2_k)
Rule('Bax4_reset',Bcl2(b=1)%Bax4(b=1) >> Bax4(b=None),deg_Bcl2_k)
Rule('MBax_reset',Bcl2(b=1)%MBax(b=1) >> MBax(b=None),deg_Bcl2_k)
bind(Bcl2,Bax4,1e-6/mito_vol,1e-3)
act3(Bax4,MOM,MOMP,1e-6/mito_vol,1e-3,1)
act2(MOMP,mCytoC,aCytoC,5e-8,1e-3,10)
act2(aCytoC,mSmac,aSmac,2e-6,1e-3,10)
migrate2(aCytoC,cCytoC,1,1e-2)
act2(cCytoC,Apaf,aApaf,5e-7,1e-3,1)
bind_consume(aApaf,C9,Apop,5e-8,1e-3)
act4(Apop,C3,C3a,5e-9,1e-3,1e-8)
migrate2(aSmac,cSmac,1,1e-2)
bind(Apop,XIAP,2e-6,1e-3)
bind(cSmac,XIAP,7e-6,1e-3)

#create_monomers(['BAD','BADp','b1433','BAD_Bcl2','BADp_Bcl2'])
#deg(BAD,3e-4)
#deg(BADp,3e-4)
#deg(b1433,3e-4)
#deg_bound(BADp,b1433,3e-4)
#deg(BAD_Bcl2,3e-4)
#Parameter('stable_BAD_Bcl2',0)
#Initial(BAD_Bcl2(b=None), stable_BAD_Bcl2)
#Observable('o_BAD_Bcl2',BAD_Bcl2(b=None))
#deg(BADp_Bcl2,3e-4)
#synth(BAD,0)
#synth(b1433,0)
#Parameter('stable_BAD',0)
#Initial(BAD(b=None),stable_BAD)
#Observable('o_BAD',BAD(b=None))
#Parameter('b1433_0',1e5)
#act2(AKTp,BAD,BADp,1e-6,1e-1,1e-1)
#act2(AKTp,BAD_Bcl2,BADp_Bcl2,1e-6,1e-1,1e-1)
#migrate(BADp,BAD,1e-3)
#bind(b1433,BADp,1e-6,1e-3)
#bind_consume(BAD,Bcl2,BAD_Bcl2,1e-6,1e-3)
#bind_consume(BADp,Bcl2,BADp_Bcl2,1e-6,1e-3)

# Necrosis
Parameter('MIM_0',1e4)
Parameter('ATP_0',1e8)
Parameter('ADP_0',5e6)
Parameter('oNa_0',1e7)
Parameter('NaPump_0',1e6)
Parameter('Mem_0',1e7)
Parameter('IP3R_0',1e5)
Parameter('eCa_0',7e6)
Parameter('PKAD_0',1e5)
Parameter('Syt7D_0',0)

# Calcium Release
create_monomers(['eCa','mCa','iCa','IP3R','IP3Ra'])
migrate2(iCa,mCa,1e-2,1e-2)
migrate(iCa,eCa,10)
act2(cCytoC,IP3R,IP3Ra,1e-6,1e-2,1)
act1(IP3Ra,eCa,iCa,1e-3)
migrate(IP3Ra,IP3R,1e-3)

# Necrosis
create_monomers(['PKAD','Syt7D','MIM','MIMP','MIMp','ROS','oNa','iNa','NaPump','Mem','dMem','Syt7a','LMem'])

act2(mCa,MIM,MIMP,1e-7,1e-2,1e-1)
migrate(MIMP,MIM,1e-2)
act2(PKAa,MIM,MIMp,2e-5,1e-1,1e-1)
migrate(MIMp,MIM,1e-3)
act1(MIM,ADP,ATP,1e4)
act1(MIMp,ADP,ATP,1e4)
synth_enz(MIMP,ROS,1e-1)
act1(MIMP,ATP,ADP,1e4)
deg(ROS,1)
migrate(ATP,ADP,4e5)
migrate(oNa,iNa,1e-2)
act2(dMem,oNa,iNa,1e-6,1e-1,10)
act4(NaPump,iNa,oNa,1e-5,1e-1,1e-8)
act3(ROS,Mem,dMem,1e-7,1e-1,1)
bind_consume(iCa,Syt,Syt7a,1e-6,1e-3)
migrate(Syt7a,iCa,1e-3)
synth_enz(Syt7a,LMem,5e-2)
deg(LMem,1e-6)
bind_consume(LMem,dMem,Mem,1e-4,0)
migrate(dMem,Mem,1e-3)


create_monomers(['TNFR2','TNFR2a'])
synth(TNFR2,0)
deg(TNFR2,1e-3)
deg(TNFR2a,1e-3)
act3(TNF,TNFR2,TNFR2a,1e-7,1e-1,1e-1)
synth(TNF,0)


######## Group Synthesis / Degradation Reactions #########
k_deg=2.9e-6

group_deg([C8,C8a,BAR,Bid,tBid,Mcl,Bax,aBax,MBax,Bax2,Bax4,MOM,MOMP,mCytoC,aCytoC,cCytoC,mSmac,aSmac,cSmac,XIAP,Apaf,aApaf,C9,Apop,C3,C3a,C3i,PARP,cPARP,C6,C6a],k_deg)

deg_PARP_k.value = 0
deg_cPARP_k.value = 0

group_degB([(FLIP,DISC),(DISC,C8),(BAR,C8a),(C8a,Bid),(Mcl,tBid),(tBid,Bax),(Bcl2,MBax),(Bcl2,Bax2),(Bcl2,Bax4),(Bax,MOM),(MOMP,mCytoC),(aCytoC,mSmac),(cSmac,XIAP),(cCytoC,Apaf),(XIAP,Apop),(Apop,C3),(C8a,C3),(XIAP,C3a),(C3a,PARP),(C3a,C6),(C6a,C8)],k_deg)

deg_C3a_PARP_k.value = 0
deg_FLIP_DISC_k.value = 1e-3

group_synth([C8,BAR,Bid,Mcl,Bax,MOM,mCytoC,mSmac,XIAP,Apaf,C9,C3,C6,PARP])


######## create initial conditions ###########
Parameter('stable_p105_p50',0)
Initial(p105_p50(s='u'),stable_p105_p50)
Observable('o_p105_p50',p105_p50(s='u'))


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

Observable('o_IKKp',IKKp(b=WILD))
Observable('o_TNF',TNF(b=WILD))
Observable('o_NFkB',p65(b=1,s='p',l='n') % p50(b=1,b_IkB=None,b_DNA=None,l='n'))
Observable('o_HSP27',HSP27n(b=WILD))
Observable('o_p38',p38p(b=WILD))
Observable('o_TNFR',TNFRa(b=WILD))

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

Observable('o_MOMP',MOMP(b=WILD))
Observable('o_iCa',iCa(b=WILD))
Observable('o_mCa',mCa(b=WILD))
Observable('o_MIMP',MIMP(b=WILD))
Observable('o_ROS',ROS(b=WILD))
Observable('o_ATP',ATP(b=WILD))
Observable('o_iNa',iNa(b=WILD))
Observable('o_cPARP',cPARP(b=WILD))
Observable('o_CytoC',cCytoC(b=WILD))
Observable('Bcl2_total',Bcl2(b=WILD))
Observable('o_IP3Ra',IP3Ra(b=WILD))
Observable('o_ADP',ADP(b=WILD))
Observable('o_tTNF',tTNF(b=WILD))
Observable('o_Syt7a',Syt7a(b=WILD))
