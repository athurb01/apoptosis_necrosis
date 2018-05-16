from pysb import *
import numpy as np

Model()

##### Initial Conditions #####

#I.C. to modulate
Parameter('PKA_0',1e5)
Parameter('Syt7_0',0)
Parameter('eCa_0',7e6)

#EARM
Parameter('L_0',3e3)
Parameter('R_0',1000)
Parameter('XIAP_0',1e5)
Parameter('Bcl2_0',3e4)
Parameter('FLIP_0',2e3)
Parameter('C8_0',1e4)
Parameter('BAR_0',1e3)
Parameter('C3_0',1e4)
Parameter('C6_0',1e4)
Parameter('PARP_0',1e6)
Parameter('mCytoC_0',5e5)
Parameter('Apaf_0',1e5)
Parameter('C9_0', 1e5)
Parameter('mSmac_0',1e5)
Parameter('Bid_0',6e4)
Parameter('Mcl_0',2e4)
Parameter('Bax_0',8e4)
Parameter('MOM_0',5e5)

#Calcium Release / Necrosis
Parameter('MIM_0',1e4)
Parameter('ATP_0',1e8)
Parameter('ADP_0',5e6)
Parameter('oNa_0',1e7)
Parameter('NaPump_0',1e6)
Parameter('Mem_0',1e7)
Parameter('IP3R_0',1e5)
Parameter('SERCA_0',1e3)

##### Compartment volume parameters #####
mito_vol = 0.01
nuc_vol = 0.3

##### Helper Functions #####

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

# 0 -> S
def synth(sub,k):
	r_name = 'synth_%s' % (sub.name)
	x = Parameter(r_name + '_k', k)
	Rule(r_name, None >> sub(b=None), x)

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

# E -> E + S
def synth_enz(enz,prod,k):
	r_name = 'synth_%s' % (prod.name)
	x = Parameter(r_name + '_k',k)
	Rule(r_name, enz() >> enz() + prod(b=None), x)

# E + S <-> E:S -> E
def deg_enz(enz,sub,kf,kr,kc):
	r1_name = '%s_%s_b' % (sub.name, enz.name)
	r2_name = 'deg_%s_by_%s' % (sub.name , enz.name)
	x = Parameter(r1_name + '_kf',kf)
	y = Parameter(r1_name + '_kr',kr)
	z = Parameter(r2_name + '_k',kc)
	Rule(r1_name, enz(b=None) + sub(b=None) | enz(b=1)%sub(b=1),x,y) 
	Rule(r2_name, enz(b=1)%sub(b=1) >> enz(b=None), z)

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


#2. Housekeeping functions

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

##### Reaction List #####
#EARM
create_monomers(['TNF','TNFR','TNFRa','XIAP','Bcl2','FLIP','BAR','Mcl','DISC','C8','C8a','C3','C3a','C3i','PARP','cPARP','C6','C6a','Bid','tBid','Bax','aBax','MBax','Bax2','Bax4','MOM','MOMP','mCytoC','aCytoC','cCytoC','mSmac','aSmac','cSmac','Apaf','aApaf','C9','Apop','ATP','ADP'])

synth(TNF,40)
deg(TNF,1e-3)
synth(TNFR,10)
deg(TNFR,1e-3)
migrate(TNFRa,DISC,1e-3)
deg(DISC,1e-3)
act3(TNF,TNFR,TNFRa,1e-7,1e-1,1e-1)
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


# Calcium Release
create_monomers(['eCa','mCa','iCa','IP3R','IP3Ra','SERCA'])
migrate2(iCa,mCa,1e-2,1e-2)
migrate(iCa,eCa,10)
act2(cCytoC,IP3R,IP3Ra,1e-6,1e-2,1)
act1(IP3Ra,eCa,iCa,1e-3)
migrate(IP3Ra,IP3R,1e-3)

# Necrosis
create_monomers(['PKA','Syt7','MIM','MIMP','MIMp','ROS','oNa','iNa','NaPump','Mem','dMem','Syt7a','LMem'])

act2(mCa,MIM,MIMP,1e-7,1e-2,1e-1)
migrate(MIMP,MIM,1e-2)
act2(PKA,MIM,MIMp,2e-5,1e-1,1e-1)
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
bind_consume(iCa,Syt7,Syt7a,1e-6,1e-3)
synth_enz(Syt7a,LMem,5e-2)
deg(LMem,1e-6)
bind_consume(LMem,dMem,Mem,1e-4,0)
migrate(dMem,Mem,1e-3)

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



##### Observables #####
Observable('o_MOMP',MOMP(b=WILD))
Observable('o_iCa',iCa(b=WILD))
Observable('o_mCa',mCa(b=WILD))
Observable('o_MIMP',MIMP(b=WILD))
Observable('o_ROS',ROS(b=WILD))
Observable('o_ATP',ATP(b=WILD))
Observable('o_iNa',iNa(b=WILD))
Observable('o_cPARP',cPARP(b=WILD))
