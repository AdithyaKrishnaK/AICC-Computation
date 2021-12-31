import pandas as pd
import math
from AICC import computeAICC

def process_data(CompPath,TempPath,VolPath):
    file = open(CompPath,'r')
    columns = file.readline().split()
    data = []
    for line in file:
        vals = line.split()
        for i in range(len(vals)):
            vals[i] = float(vals[i])

        containment = parse(columns,vals)
        add_temp(containment,TempPath,vals[0])
        add_vol(containment,VolPath)
        if not validate(containment):
            print("Please verify the data")
            return
        
        res = {"TIME":vals[0]}
        for space in containment:
            [T,P] = computeAICC(space)
            res[space.get('ID')+'_T'] = T
            res[space.get('ID')+'_P'] = P
        
        data.append(res)
    
    result = pd.DataFrame.from_dict(data)
    result.to_csv('result.csv')
        



def search_by_id(listName,spaceID):
    for i in listName:
        if i['ID'] == spaceID:
            return True
    
    return False

def add_param(listName,spaceID,param,val):
    for i in listName:
        if i['ID'] == spaceID:
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
            [comp,space] = parameter.split(',')
            if search_by_id(containment,space):
                add_param(containment,space,comp[2:],val)
            else:
                create_space(containment, space)
                add_param(containment,space,comp[2:],val)

        elif parameter.startswith('P_TOT'):
            space = parameter.replace('P_TOT','')
            if search_by_id(containment,space):
                add_param(containment,space,"P",val)
            else:
                create_space(containment, space)
                add_param(containment,space,"P",val)
    
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
                        add_param(containment,col[2:],"T",vals[i])

def add_vol(containment, path,time):
    file = open(path,'r')
    columns = file.readline().split()
    vals = file.readline().split()
    for i in range(len(columns)):
        if columns[i].startswith('V'):
            col = columns[i]
            if search_by_id(containment,col[2:]):
                add_param(containment,col[2:],"V",vals[i])


def validate(containment):
    flg = True
    tol = math.exp(10,-5)
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
        
        if total_x>1:
            print("Sum of mole fractions cannot be greater than 1 in", space.get('ID'))
            flg = False
        elif 1-total_x < tol:
            print("Sum of mole fractions have to be one. Data missing in", space.get('ID'))
            flg = False

    return flg

