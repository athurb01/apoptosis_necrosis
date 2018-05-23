### Written by Krista Pullen 
### Last Updated 12/20/2017
###Figure 7(a). Requires model_Syt7.py and run_model_Syt7.py

from run_model_Syt7 import *
import matplotlib.pyplot as plt

chg_val(m.PKA_0,0) #figure 7a compares the relationship between mitochondrial calcium loading and MPT by scanning abundances of ER calcium stores, so I just turn off any inhibitory effect of PKA here
W,X,Y,Z = create_phase() #you'll notice this function is really just designed for this figure and 7c
pl.ion()
fig = pl.figure()
pl.plot(W,X)
plt.xlabel('log Mito Ca++')
plt.ylabel('% MPT')

###Figure 7(b). Requires model_Syt7.py and run_model_Syt7.py

from run_model_Syt7 import *

chg_val(m.eCa_0,7e6) #if you've just done 7a and scanned ER calcium store abundances up to 1e8, this resets to default value

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X)

plt.xlabel('log PKA (Copies)')

plt.ylabel('% MPT')

###Figure 7(c). Requires model_Syt7.py and run_model_Syt7.py

from run_model_Syt7 import *

chg_val(m.synth_ROS_k,0.1)

W,X,Y,Z = create_phase()

pl.ion()

pl.figure()

pl.plot(X,Z)

chg_val(m.synth_ROS_k,0)

W,X,Y,Z = create_phase()

pl.plot(X,Z)

chg_val(m.synth_ROS_k,0.05)

W,X,Y,Z = create_phase()

pl.plot(X,Z)

chg_val(m.synth_ROS_k,0.01)

W,X,Y,Z = create_phase()

pl.plot(X,Z)

plt.xlabel('% MPT')

plt.ylabel('Osmotic Stress (fold [Na] increase from baseline)')

plt.legend(['Mem Damage (ROS_ks = 0.1)','No Mem Damage (ROS_ks = 0)','ROS_ks = 0.05','ROS_ks = 0.01'])

###Figure 7(d). Requires model_Necrosis.py and run_model_Necrosis.py

from run_model_Necrosis import *

chg_val(m.Syt7_0,0)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.figure()

pl.plot(a,Z)

chg_val(m.Syt7_0,2500)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,Z)

chg_val (m.Syt7_0,5000)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,Z)

chg_val (m.Syt7_0,7500)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,Z)

chg_val(m.Syt7_0,10000)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,Z)

plt.scatter(a,Z,s=None,c='k',marker='X')

plt.legend(['Syt7 = 0','Syt7 = 2500','Syt7 =5000','Syt7 =7500','Syt7 =10000','No Mem Damage'])

plt.xlabel('Log PKA (Copies)')

plt.ylabel('Osmotic Stress (Fold [Na] Increase from Baseline)')

###Figure 8. Requires model_Necrosis.py and run_model_Necrosis.py

from run_model_Necrosis import *

pl.figure()

fill_between([0,2], 0, 100, where=None, interpolate=True, step=None, data=None, color = 'r')

fill_between([2,4], 0, 100, where=None, interpolate=True, step=None, data=None, color = 'y')

fill_between([4,5], 0, 100, where=None, interpolate=True, step=None, data=None, color = 'g')

plt.xlim([0,5])

plt.ylim([0,100])

chg_val(m.Syt7_0,0)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

plt.scatter(a,Z,s=None,c='k',marker='+')

chg_val(m.Syt7_0,10000)

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

plt.scatter(a,Z,s=None,c='k',marker='x')

plt.axhline(y=40, color='k', linestyle='-')

plt.axvline(x=2, color='k', linestyle='-')

plt.axvline(x=4, color='k', linestyle='-')

###Figure 9(a). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

Chg_val(m.LO_AA_b_kf, 1e-08)

scan_custom(m.p38p_HSP27_b_kf,[-7,-5],'o_Cox')

###Figure 9(b). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(p38p_HSP27_b_kf,3e-6)

scan_custom(m.synth_PIP3_k,[-4,0],'o_Cox')

###Figure 10(a). Requires Frontside.py and run_Frontside.py

from run_Frontside import *

chg_val(m.LO_AA_b_kf,1e-5)

chg_val(m.LO_0,1e5)

graph('o_LXA4')

chg_val(m.LO_0,1e4)

graph('o_LXA4',prep_fig=False)

chg_val(m.LO_0,5e3)

graph('o_LXA4',prep_fig=False)

chg_val(m.LO_0,1e3)

###Figure 10(b). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_0,1e3)

pl.figure()

s.run()

output = s.yobs['o_PGE2'] + s.yobs['o_EP2'] + s.yobs['o_EP4']

pl.plot(t/3600,output)

chg_val(m.LO_0,5e3)

s.run()

output = s.yobs['o_PGE2'] + s.yobs['o_EP2'] + s.yobs['o_EP4']

pl.plot(t/3600,output)

chg_val(m.LO_0,1e4)

s.run()

output = s.yobs['o_PGE2'] + s.yobs['o_EP2'] + s.yobs['o_EP4']

pl.plot(t/3600,output)

chg_val(m.LO_0,1e5)

s.run()

output = s.yobs['o_PGE2'] + s.yobs['o_EP2'] + s.yobs['o_EP4']

pl.plot(t/3600,output)
graph('o_LXA4',prep_fig=False)

###Figure 10(c). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_0,1e3)

graph('o_NFkB')

chg_val(m.LO_0,5e3)

graph('o_NFkB',prep_fig=False)

chg_val(m.LO_0,1e4)

graph('o_NFkB',prep_fig=False)

chg_val(m.LO_0,1e5)

graph('o_NFkB',prep_fig=False)

###Figure 11. Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

LO_wrap([3,5],'o_NFkB',norm=True)

comp_LO_TNF([3,5],norm=True)

LO_wrap([3,5],'o_Cox',norm=True)

LO_wrap([3,5],'o_Syt7',norm=True)

comp_LO_PGE2([3,5],norm=True)

LO_wrap([3,5],'o_PKA',norm=True)

LO_wrap([3,5],'o_DISC',norm=True)

###Figure 12(a). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.synth_AA_k, 0.05)

comp_LO_PGE2([3,5])

chg_val(m.synth_AA_k, 0.01)

comp_LO_PGE2([3,5])

chg_val(m.synth_AA_k, 0.002)

comp_LO_PGE2([3,5])

###Figure 12(b). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.synth_AA_k, 0.05)

LO_wrap([3,5],'o_LXA4')

chg_val(m.synth_AA_k, 0.01)

LO_wrap([3,5],'o_LXA4')

chg_val(m.synth_AA_k, 0.002)

LO_wrap([3,5],'o_LXA4')

###Figure 12(c). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.synth_AA_k, 0.05)

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_AA_k, 0.01)

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_AA_k, 0.002)

LO_wrap([3,5],'o_PKA')

###Figure 12(d). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.synth_AA_k, 0.05)

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_AA_k, 0.01)

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_AA_k, 0.002)

LO_wrap([3,5],'o_Syt7')

###Figure 13. Requires FullModel.py and run_FullModel.py

from run_FullModel.py import *

Osm_LO([3,5])

###Figure 14. Requires FullModel.py and run_FullModel.py

from run_FullModel.py import *

chg_val(m.Pol_dBcl2_p_kf,1e-8)

time_to_nec([3.0,5.0],'b')

chg_val(m.Pol_dBcl2_p_kf,3e-8)

time_to_nec([3.0,5.0],'g')

chg_val(m.Pol_dBcl2_p_kf,5e-8)

time_to_nec([3.0,5.0],'r')

###Figure 15(a). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_AA_b_kf,1e-5)

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_TNF_k,60)

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_TNF_k,0)

chg_val(m.synth_TNFR2_k, 60)

LO_wrap([3,5],'o_PKA')

###Figure 15(b). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_AA_b_kf,1e-5)

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_TNF_k,60)

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_TNF_k,0)

chg_val(m.synth_TNFR2_k, 60)

LO_wrap([3,5],'o_Syt7')

###Figure 15(c). Requires FullModel.py and run_FullModel.py

from run_FullModel.py import *

chg_val(m.LO_AA_b_kf,1e-5)

time_to_nec([3.0,5.0],'b')

chg_val(m.synth_TNF_k,60)

time_to_nec([3.0,5.0],'r')

chg_val(m.synth_TNF_k,0)

chg_val(m.synth_TNFR2_k, 60)

time_to_nec([3.0,5.0],'g')

###Figure 16(a). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_AA_b_kf,1e-5)

chg_val(m.synth_p50_k,45)

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_p50_k,75)

stable_kB()

LO_wrap([3,5],'o_PKA')

chg_val(m.synth_p50_k,30)

stable_kB()

LO_wrap([3,5], 'o_PKA')

###Figure 16(b). Requires Frontside.py and run_Frontside.py

from run_Frontside.py import *

chg_val(m.LO_AA_b_kf,1e-5)

chg_val(m.synth_p50_k,45)

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_p50_k,75)

stable_kB()

LO_wrap([3,5],'o_Syt7')

chg_val(m.synth_p50_k,30)

stable_kB()

LO_wrap([3,5],'o_Syt7')

###Figure 16(c). Requires FullModel.py and run_FullModel.py

from run_FullModel.py import *

chg_val(m.LO_AA_b_kf,1e-5)

chg_val(m.synth_p50_k,45)

time_to_nec([3.0,5.0],'b')

chg_val(m.synth_p50_k,75)

stable_kB()

time_to_nec([3.0,5.0],'g')

chg_val(m.synth_p50_k,30)

stable_kB()

time_to_nec([3.0,5.0],'r')

##############################################################################################


###Figure 17. Requires model_Syt7.py and run_model_Syt7.py. 

from run_model_Syt7 import *

chg_val(m.eCa_0,1e3) 

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X,'k')

chg_val(m.eCa_0,1e4) 

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X,'b')

chg_val(m.eCa_0,1e5) 

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X,'g')

chg_val(m.eCa_0,1e6) 

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X,'r')

chg_val(m.eCa_0,1e7) 

chg_val(m.PKA_MIM_b_kf,2e-5) 

a,X,Y,Z = create_phase2(m.PKA_0,[0,5])

pl.plot(a,X,'y')

plt.xlabel('log PKA (Copies)')

plt.ylabel('% MPT')

###Figure 17. Requires FullModel.py and run_FullModel.py

from run_FullModel.py import *

chg_val(m.Pol_dBcl2_p_kf,3e-8)
time_to_nec([3.0,6.0],'g', 1000000.0, 1000000.0)
time_to_nec([3.0,6.0],'k', 750000.0, 1000000.0)
time_to_nec([3.0,6.0], 'b', 500000.0, 1000000.0)
time_to_nec([3.0,6.0],'k', 250000.0, 1000000.0)