import Frontside as m
import pylab as pl
import numpy as np
from pysb.integrate import *
tmax = 24
t = pl.linspace(0, 3600 * tmax)
s = Solver(m.model, t)
s.run()

def qparams(list):
    for i in m.model.parameters:
        tag = True
        for j in list:
            if j not in i.name:
                tag = False

        if tag:
            print i


def chg_val(param, k):
    param.value = k


def graph(obs, prep_fig = True):
    s.run()
    if prep_fig:
        pl.ion()
        pl.figure()
    pl.plot(t / 3600, s.yobs[obs])


def stable_kB(init_IKK = 100000.0, run = False):
    chg_val(m.TB_0, 0)
    list = []
    list2 = []
    for i in m.model.parameters:
        if 'stable_' in i.name:
            chg_val(i, 0)
            list.append(i)
            list2.append(i.name[7:])

    chg_val(m.stable_IKK, init_IKK)
    chg_val(m.TB_0, 0)
    chg_val(m.AKTp_IKK_b_kf, 0)
    chg_val(m.TNFRa_TAK_b_kf, 0)
    s.run()
    for i in range(len(list)):
        chg_val(list[i], s.yobs['o_' + list2[i]][-1])

    chg_val(m.TB_0, 1)
    chg_val(m.AKTp_IKK_b_kf, 1e-05)
    chg_val(m.TNFRa_TAK_b_kf, 1e-07)
    if run:
        s.run()

def scan_custom(param,range,obs,norm=False):
    if norm:
                chg_val(param,0)
                s.run()
                Pmax = s.yobs[obs][-1]
    else:
                Pmax = 1
    obs_out = []
    logrange = pl.linspace(range[0],range[1])
    linrange = 10**logrange
    for r in linrange:
                chg_val(param,r)
                s.run()
                obs_out.append(s.yobs[obs][-1]/Pmax)
    pl.plot(logrange,obs_out,label=obs[2:])
        

def scan_LXA4(range, obs, norm = False):
    if norm:
        chg_val(m.LO_AA_b_kf, 0)
        s.run()
        Pmax = s.yobs[obs][-1]
    else:
        Pmax = 1
    obs_out = []
    for r in range:
        chg_val(m.LO_AA_b_kf, r)
        s.run()
        obs_out.append(s.yobs[obs][-1] / Pmax)

    return obs_out


def composite_TNF(range, norm = False):
    if norm:
        chg_val(m.LO_AA_b_kf, 0)
        s.run()
        Pmax = s.yobs['o_TNF'][-1] + s.yobs['o_TNFR'][-1]
    else:
        Pmax = 1
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    TNF_out = scan_LXA4(linrange, 'o_TNF')
    TNFR_out = scan_LXA4(linrange, 'o_TNFR')
    TNF_array = np.array(TNF_out)
    TNFR_array = np.array(TNFR_out)
    Comp_array = (TNF_array + TNFR_array) / Pmax
    pl.plot(logrange, Comp_array, label='TNF')


def composite_PGE2(range, norm = False):
    if norm:
        chg_val(m.LO_AA_b_kf, 0)
        s.run()
        Pmax = s.yobs['o_PGE2'][-1] + s.yobs['o_EP2'][-1] + s.yobs['o_EP4'][-1]
    else:
        Pmax = 1
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    PGE2_out = scan_LXA4(linrange, 'o_PGE2')
    EP2_out = scan_LXA4(linrange, 'o_EP2')
    EP4_out = scan_LXA4(linrange, 'o_EP4')
    Comp_array = (np.array(PGE2_out) + np.array(EP2_out) + np.array(EP4_out)) / Pmax
    pl.plot(logrange, Comp_array, label='PGE2')


def scan_wrap(range, obs, norm = False):
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    obs_out = scan_LXA4(linrange, obs, norm)
    LXA4_out = scan_LXA4(linrange, 'o_LXA4')
    pl.plot(LXA4_out, obs_out, label=obs[2:])


def PGE2_LXA4(range):
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    PGE2_out = scan_LO(linrange, 'o_PGE2')
    EP2_out = scan_LO(linrange, 'o_EP2')
    EP4_out = scan_LO(linrange, 'o_EP4')
    LXA4_out = scan_LO(linrange, 'o_LXA4')
    Comp_array = np.array(PGE2_out) + np.array(EP2_out) + np.array(EP4_out)
    LXA4_array = np.array(LXA4_out)
    ratio = np.log10(Comp_array / LXA4_array)
    pl.plot(logrange, ratio, label='PGE2/LXA4 ratio')


def scan_LO(range, obs, norm = False):
    if norm:
        chg_val(m.LO_0, 0)
        s.run()
        Pmax = s.yobs[obs][-1]
    else:
        Pmax = 1
    obs_out = []
    for r in range:
        chg_val(m.LO_0, r)
        s.run()
        obs_out.append(s.yobs[obs][-1] / Pmax)

    return obs_out


def LO_wrap(range, obs, norm = False):
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    obs_out = scan_LO(linrange, obs, norm)
    pl.plot(logrange, obs_out, label=obs[2:])


def comp_LO_TNF(range, norm = False):
    if norm:
        chg_val(m.LO_0, 0)
        s.run()
        Pmax = s.yobs['o_TNF'][-1] + s.yobs['o_TNFR'][-1]
    else:
        Pmax = 1
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    TNF_out = scan_LO(linrange, 'o_TNF')
    TNFR_out = scan_LO(linrange, 'o_TNFR')
    TNF_array = np.array(TNF_out)
    TNFR_array = np.array(TNFR_out)
    Comp_array = (TNF_array + TNFR_array) / Pmax
    pl.plot(logrange, Comp_array, label='TNF')


def comp_LO_PGE2(range, norm = False):
    if norm:
        chg_val(m.LO_0, 0)
        s.run()
        Pmax = s.yobs['o_PGE2'][-1] + s.yobs['o_EP2'][-1] + s.yobs['o_EP4'][-1]
    else:
        Pmax = 1
    logrange = pl.linspace(range[0], range[1])
    linrange = 10 ** logrange
    PGE2_out = scan_LO(linrange, 'o_PGE2')
    EP2_out = scan_LO(linrange, 'o_EP2')
    EP4_out = scan_LO(linrange, 'o_EP4')
    Comp_array = (np.array(PGE2_out) + np.array(EP2_out) + np.array(EP4_out)) / Pmax
    pl.plot(logrange, Comp_array, label='PGE2')


stable_kB()