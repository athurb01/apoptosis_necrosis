# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 20:41:28 2018

@author: Amy Thurber
"""

from pysb import *
from pysb.macros import *
import numpy as np

model = Model('pysb_explore')

Monomer('A', ['b'])
Monomer('B', ['b'])

Observable('o_A', A(b=None))
Observable('o_B', B(b=None))
Observable('o_AB', A(b=1)%B(b=1))

Parameter('A_0', 100)
Parameter('B_0', 150)

Parameter('binds', 0.5)
Parameter('degrades', 0.1)

Initial(A(b=None), A_0)
Initial(B(b=None), B_0)

Rule('binding', A(b=None) + B(b=None) >> A(b=1)%B(b=1), binds)
Rule('degrading', A(b=ANY) >> None, degrades)
