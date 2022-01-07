from functools import partial
import pandas as pd
import time
from AICC import computeAICC, findEqT
import pyromat as pm

not_composition = ['P','V','T','ID']
ignore_compartments = ['SCONTAIN','ENVIRON']

def computeU(r,n_r,t):
  pm.config['unit_matter'] = 'kmol'
  data = {}
  R = 8.314
  for key in r.keys():
    if key not in not_composition:
      data[key] = pm.get('ig.'+key)
  val = 0
  for key in r.keys():
    if key not in not_composition:
      val += data[key].e(T=t)*r[key]*n_r
  return val

def homogenise(containment):
    total = {}
    enthalpy = 0
    for space in containment:
            if space.get('ID').upper() not in ignore_compartments:
                for key in space.keys():
                    if key not in not_composition:
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

    fun = partial(energy_balance,total,n_tot,enthalpy)
    tf = findEqT(fun,315)

    total['T'] = tf
    total['P'] = n_tot*total['T']*8.314/total['V']
    return computeAICC(total)

def energy_balance(total,n,enthalpy,t):
    return computeU(total,n,t) - enthalpy

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
            if key not in not_composition:
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
                if key not in not_composition:
                    tot += space[key]
            
            for key in space.keys():
                if key not in not_composition:
                    space[key] /= tot
        
        if not validate(containment):
            print(containment)
            print("Please verify the data")
            return
        
        b =time.time()
        #print("Parsing",b-a)
        res = {"TIME":vals[0]}
        for space in containment:
            if space.get('ID').upper() not in ignore_compartments:
                ans = computeAICC(space)
                res[space.get('ID')+'_T'] = float(ans['T'])
                res[space.get('ID')+'_P'] = float(ans['P'])
        
        
        ans = homogenise(containment)
        res['Homogenised'+'_T'] = float(ans['T'])
        res['Homogenised'+'_P'] = float(ans['P'])
        data.append(res)
        
        a =time.time()
        #print("Proccessing",a-b)

    result = pd.DataFrame.from_dict(data)
    result.to_csv('result.csv')

# begin = time.time()
# process_data('data/H2-2.plot','data/temperature-flu.plot','data/vol.txt')
# end = time.time()
# print("Time taken is ", end - begin)