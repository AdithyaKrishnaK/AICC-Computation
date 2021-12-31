import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import math
import pyromat as pm

def uchange(r,p,T0,t):
    pm.config['unit_matter'] = 'kmol'
    o2 = pm.get('ig.O2')
    n2 = pm.get('ig.N2')
    h2 = pm.get('ig.H2')
    h2o = pm.get('ig.H2O')

    val = p['H2O']*h2o.e(T=t) + p['O2']*o2.e(T=t) + p['N2']*n2.e(T=t) + p['H2']*h2.e(T=t) - r['H2']*h2.e(T=T0) - r['O2']*o2.e(T=T0) - r['N2']*n2.e(T=T0) 
    return val


def findEqT(r,p,T0):
    a = T0
    b = 10000
    while abs((b-a)/b)>0.00000001:
        c = (b+a)/2
        if uchange(r,p,T0,c)<0:
            a = c
        else:
            b = c
    
    return a

def computeAICC(r):
    P0 = r['P']
    T0 = r['T']
    V = r['V']
    p = {}
    if r['H2']/2>r['O2']:
        p['H2O'] = r['O2']*2
        p['H2']= r['H2']-r['O2']*2
        p['O2'] = 0
    else:
        p['H2O'] = r['H2']
        p['H2']= 0
        p['O2'] = r['O2'] - r['H2']/2

    p['n2'] = r['n2']
    n_p = p['H2']+p['O2']+p['H2O']+p['N2']

    R = 8.314

    T = findEqT(r,p,T0)
    P_AICC = n_p*R*T/V
    return [T,P_AICC]