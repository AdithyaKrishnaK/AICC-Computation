import math
import numpy as np
import matplotlib.pyplot as plt


def calcVf(x_dry,spray,x_d,L):
    vf = 0
    if x_dry<=0.1:
        a = -4.877
        b = -3.008
        if spray:
            vf = math.pow(L,0.333)*(59.65*x_dry-1.248)*math.exp(x_d*(a+b*x_d))
        else:
            vf = math.pow(L,0.333)*(23.70*x_dry-0.862)*math.exp(x_d*(a+b*x_d))
    elif x_dry<=0.18:
        a = -4.877
        b = -3.008
        if spray:
            vf = math.pow(L,0.333)*(2074*x_dry*x_dry-347.23*x_dry+18.7)*math.exp(x_d*(a+b*x_d))
        else:
            vf = math.pow(L,0.333)*(1724*x_dry*x_dry-267.28*x_dry+10.996)*math.exp(x_d*(a+b*x_d))
    elif x_dry<0.25:
        vf = (x_dry*(calcVf(0.18,spray,x_d,L)- calcVf(0.25,spray,x_d,L)) - (calcVf(0.18,spray,x_d,L)*0.25- calcVf(0.25,spray,x_d,L)*0.18))/(0.18-0.25)
    elif x_dry<=0.35:
        a = -0.641
        b = -18.38
        vf = math.pow(L,0.333)*(289.73*x_dry-33.769)*math.exp(x_d*(a+b*x_d))
    elif x_dry<0.45:
        vf = (x_dry*(calcVf(0.45,spray,x_d,L)- calcVf(0.35,spray,x_d,L)) - (calcVf(0.45,spray,x_d,L)*0.35- calcVf(0.35,spray,x_d,L)*0.45))/(0.45-0.35)
    elif x_dry<=0.72:
        a = -17.279
        b = 18.07
        vf = math.pow(L,0.333)*(-199.62*x_dry+145.07)*math.exp(x_d*(a+b*x_d))
    
    return vf

def flame(r, vol, spray, data = {'F':0.541,'x_c':0.07,'x_o2':0.05,'x_d':0.55}):
    R = 8.314
    n_r = r['P']*r['V']/(R*r['T'])

    if r.get('H2') != None:
        h_key = 'H2'
    else:
        h_key = 'D2' 

    if r.get('CO') == None:
        r['CO'] = 0

    x_c = r[h_key] + data['F']*r['CO']
    
    x_d = 0
    for key in r.keys():
        if key in ['CO2','D2O','H2O']:
            x_d += r[key]

    if r.get('O2')==None:
        print("No oxygen. No explosion")
        return

    deltaN = 0
    if r.get('N2') != None:
        deltaN = 0.79*(r['N2']-3.774*r['O2'])
        x_d += max(deltaN,0)

    if x_c>=data['x_c'] and r['O2']>=data['x_o2'] and x_d<=data['x_d']:
        print("No Explosion")
        return
    L =  math.pow(vol,0.3333)
    if r.get('H2O') == None:
        r['H2O'] = 0
    if r.get('CO2') == None:
        r['CO2'] = 0
    x_dry = x_c/(1-r['H2O']-r['CO2']-max(deltaN,0))
    vf = calcVf(x_dry,spray,x_d,L)


x = np.linspace(0.1,0.72,num=62)
vf = []
for i in x:
    o =  (1-i)/4.76
    x_d = 3.76*o
    vf.append(calcVf(i,False,0,1))

plt.plot(x,vf)
plt.show()