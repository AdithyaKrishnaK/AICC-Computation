from typing import Union
import pyromat as pm

not_composition = ['P','V','T','ID']

def uchange(r,p,t,n=0):
  pm.config['unit_matter'] = 'kmol'
  data = {}
  R = 8.314
  if n==0:
    n_r = r['P']*r['V']/(R*r['T'])
  else:
    n_r = n
  for key in r.keys():
    if key not in not_composition:
      data[key] = pm.get('ig.'+key)
  for key in p.keys():
    if key not in not_composition and data.get(key) == None:
      data[key] = pm.get('ig.'+key)
  val = 0
  for key in r.keys():
    if key not in not_composition:
      val -= data[key].e(T=r['T'])*r[key]*n_r
  
  for key in p.keys():
    if key not in not_composition:
      val += data[key].e(T=t)*p[key]
  return val


def findEqT(r,p):
  a = r['T']
  b = 4000
  while abs((b-a)/b)>0.00000001:
    c = (b+a)/2
    if uchange(r,p,c)*uchange(r,p,b)<0:
      a = c
    else:
      b = c
  
  return a



def computeProd(r,data):
  p = {}
  R = 8.314
  n_r = r['P']*r['V']/(R*r['T'])

  if r.get('H2') != None:
    h_key = 'H2'
  else:
    h_key = 'D2' 

  if r.get('CO') == None:
    r['CO'] = 0

  x = r[h_key] + data['F']*r['CO']
  
  x_d = 0
  for key in r.keys():
    if key in ['CO2','D2O','H2O']:
      x_d += r[key]

  if r.get('N2') != None and r.get('O2')!= None:
      x_d += 0.79*(r['N2']-3.774*r['O2'])

  if x>=data['x_c'] and r['O2']>=data['x_o2'] and x_d<=data['x_d']:
    if (r[h_key]+r['CO'])/2>r['O2']:
        k = 2*r['O2']/(r[h_key]+r['CO']) 
        p[h_key+'O'] = k*r[h_key]*n_r
        p['CO2'] = k*r['CO']*n_r
        p['CO'] = (r['CO']-k*r['CO'])*n_r
        p[h_key]= (r[h_key]-k*r[h_key])*n_r
        if p[h_key] < 0 or p['CO']<0:
          print("Error")
        p['O2'] = 0
    else:
        p[h_key+'O'] = r[h_key]*n_r
        p[h_key]= 0
        p['CO'] = 0
        p['CO2'] = r['CO']*n_r
        p['O2'] = (r['O2'] - (r[h_key]+r['CO'])/2)*n_r
  else:
    for key in r.keys():
        if key not in not_composition:
            p[key] = n_r*r[key]
    return p
  
  for key in r.keys():
    if key not in Union(['H2','D2','O2','D2O','H2O','CO','CO2'],not_composition):
      p[key] = n_r*r[key]
  return p

def computeAICC(r, data = {'F':0.541,'x_c':0.07,'x_o2':0.05,'x_d':0.55}):
  
  if r.get('H2') == None and r.get('D2') == None:
    print("No hydrogen or deutrium")
    return

  R = 8.314
  p = computeProd(r,data)

  n_p=0
  for key in p.keys():
    if key not in not_composition:
      n_p += p[key]

  T = findEqT(r,p)
  P_AICC = n_p*R*T/r['V']
  p['T'] = T
  p['P'] = P_AICC
  return p