import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from AICC import computeAICC, uchange
import pyromat as pm

def computeU(r,n_r,t,enthalpy = 0):
  pm.config['unit_matter'] = 'kmol'
  data = {}
  R = 8.314
  for key in r.keys():
    if key not in ['P','T','V','ID']:
      data[key] = pm.get('ig.'+key)
  val = 0
  for key in r.keys():
    if key not in ['P','T','V','ID']:
      val += data[key].e(T=t)*r[key]*n_r
  return val - enthalpy

def homogenise(containment):
    total = {}
    enthalpy = 0
    for space in containment:
            if space.get('ID').upper() not in ['SCONTAIN','ENVIRON']:
                for key in space.keys():
                    if key not in ['P','T','V','ID']:
                        if total.get(key) != None:
                            total[key] = total[key] + space[key]*space['P']*space['V']/(8.314*space['T'])
                        else:
                            total[key] = space[key]*space['P']*space['V']/(8.314*space['T'])
                    if key == 'V':
                        if total.get(key) != None:
                            total[key] = total[key] + space[key]
                        else:
                            total[key] = space[key]
                enthalpy += computeU(space,space['P']*space['V']/(8.314*space['T']),space['T'])
    n_tot = 0
    for key in total:
        if key != 'V':
            n_tot += total[key]
    for key in total:
        if key != 'V':
            total[key] /= n_tot

    a = 0
    b = 1000
    while abs((b-a)/b)>0.00000001:
        c = (b+a)/2
        if float(computeU(total,n_tot,c,float(enthalpy)))*float(computeU(total,n_tot,b,float(enthalpy)))<0:
            a = c
        else:
            b = c
    
    # a = np.linspace(100,10000,num=4000)
    
    # plt.plot(a,computeU(total,n_tot,a,float(enthalpy)))
    # plt.show()
    if float(computeU(total,n_tot,a,float(enthalpy))) > 0.000001:
        print("Not converged")
        exit()
    total['T'] = a
    total['P'] = n_tot*total['T']*8.314/total['V']
    return computeAICC(total)

def search_by_id(listName,spaceID):
    for i in listName:
        if i['ID'].upper() == spaceID.upper():
            return True
    
    return False

def add_param(listName,spaceID,param,val):
    for i in listName:
        if i['ID'].upper() == spaceID.upper():
            i[param] = val
            
def create_space(listName,spaceID):
    new_obj = {'ID':spaceID}
    listName.append(new_obj)

def parse(col,row):
    containment = []
    for i in range(len(col)):
        parameter = col[i]
        val = row[i]
        if parameter.startswith('PP'):
            [comp,space] = parameter.split('_')
            if search_by_id(containment,space):
                add_param(containment,space,comp[2:],float(val))
            else:
                create_space(containment, space)
                add_param(containment,space,comp[2:],float(val))

        elif parameter.startswith('P_TOT'):
            space = parameter.replace('P_TOT','')
            space = space[1:] if space.startswith('_') else space
            if search_by_id(containment,space):
                add_param(containment,space,"P",float(val))
            else:
                create_space(containment, space)
                add_param(containment,space,"P",float(val))
    
    return containment

def add_temp(containment, path,time):
    file = open(path,'r')
    columns = file.readline().split()
    for line in file:
        vals = line.split()
        if float(vals[0]) == time:
            for i in range(len(columns)):
                if columns[i].startswith('T'):
                    col = columns[i]
                    if search_by_id(containment,col[2:]):
                        add_param(containment,col[2:],"T",float(vals[i]))

def add_vol(containment, path):
    file = open(path,'r')
    for line in file:
        vals = line.split()
        if vals[0].startswith('V'):
            col = vals[0]
            if search_by_id(containment,col[2:]):
                add_param(containment,col[2:],"V",float(vals[1]))


def validate(containment):
    flg = True
    tol = 0.000001
    for space in containment:
        if space.get('T') == None:
            print("Temperature is missing in ", space.get("ID"))
            flg = False
        if space.get('P') == None:
            print("Pressure is missing in ", space.get("ID"))
            flg = False
        if space.get('V') == None:
            print("Volume is missing in ", space.get("ID"))
            flg = False
        total_x = 0
        for key in space.keys():
            if key not in ['ID','T','P','V']:
                total_x += space.get(key)
        
        if total_x>1+tol:
            print("Sum of mole fractions cannot be greater than 1 in", space.get('ID'))
            print(total_x)
            flg = False
        elif 1-total_x > tol:
            print("Sum of mole fractions have to be one. Data missing in", space.get('ID'))
            flg = False
            print(1-total_x)

    return flg


def process_data(CompPath,TempPath,VolPath):
    file = open(CompPath,'r')
    columns = file.readline().split()
    data = []
    for line in file:
        a= time.time()
        vals = line.split()
        for i in range(len(vals)):
            vals[i] = float(vals[i])
        containment = parse(columns,vals)
        add_temp(containment,TempPath,vals[0])
        add_vol(containment,VolPath)
        
        for space in containment:
            tot = 0
            for key in space.keys():
                if key not in ['ID','V','P','T']:
                    tot += space[key]
            
            for key in space.keys():
                if key not in ['ID','V','P','T']:
                    space[key] /= tot
        
        if not validate(containment):
            print(containment)
            print("Please verify the data")
            return
        
        b =time.time()
        print("Parsing",b-a)
        res = {"TIME":vals[0]}
        for space in containment:
            if space.get('ID').upper() not in ['SCONTAIN','ENVIRON']:
                ans = computeAICC(space)
                res[space.get('ID')+'_T'] = ans['T']
                res[space.get('ID')+'_P'] = ans['P']
        
        
        ans = homogenise(containment)
        res['Homogenised'+'_T'] = ans['T']
        res['Homogenised'+'_P'] = ans['P']
        data.append(res)
        
        a =time.time()
        print("Proccessing",a-b)

    result = pd.DataFrame.from_dict(data)
    result.to_csv('result.csv')

# begin = time.time()
# process_data('data/H2-2.plot','data/temperature-flu.plot','data/vol.txt')
# end = time.time()
# print("Time taken is ", end - begin)