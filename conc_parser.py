from functools import partial
import pandas as pd
from AICC import computeAICC, inbuiltsolver, not_composition
import pyromat as pm
from tqdm import tqdm

ignore_compartments = ["SCONTAIN", "ENVIRON"]
column_prefix = {"partialP": "PP", "totalP": "P_TOT"}


def computeU(r, n_r, t):
    pm.config["unit_matter"] = "kmol"
    data = {}
    R = 8.314
    for key in r.keys():
        if key not in not_composition:
            data[key] = pm.get("ig." + key)
    val = 0
    for key in r.keys():
        if key not in not_composition:
            val += data[key].e(T=t) * r[key] * n_r
    return val


def homogenise(containment):
    total = {}
    enthalpy = 0
    maxT = 0
    minT = 6000
    for space in containment:
        if space.get("ID").upper() not in ignore_compartments:
            if space.get("T") > maxT:
                maxT = space.get("T")
            if space.get("T") < minT:
                minT = space.get("T")
            for key in space.keys():
                if key not in not_composition:
                    if total.get(key) != None:
                        total[key] = total[key] + space[key] * space["P"] * space[
                            "V"
                        ] / (8.314 * space["T"])
                    else:
                        total[key] = (
                            space[key] * space["P"] * space["V"] / (8.314 * space["T"])
                        )
                if key == "V":
                    if total.get(key) != None:
                        total[key] = total[key] + space[key]
                    else:
                        total[key] = space[key]
            enthalpy += computeU(
                space, space["P"] * space["V"] / (8.314 * space["T"]), space["T"]
            )
    n_tot = 0
    for key in total:
        if key != "V":
            n_tot += total[key]
    for key in total:
        if key != "V":
            total[key] /= n_tot

    fun = partial(energy_balance, total, n_tot, enthalpy)
    tf = inbuiltsolver(fun, minT, maxT)

    total["T"] = tf
    total["P"] = n_tot * total["T"] * 8.314 / total["V"]
    return computeAICC(total)


def energy_balance(r, n, enthalpy, t):
    return computeU(r, n, t) - enthalpy


def search_by_id(listName, spaceID):
    for i in listName:
        if i["ID"].upper() == spaceID.upper():
            return True

    return False


def add_param(listName, spaceID, param, val):
    for i in listName:
        if i["ID"].upper() == spaceID.upper():
            i[param] = val


def create_space(listName, spaceID):
    new_obj = {"ID": spaceID}
    listName.append(new_obj)


def parse(containment, col, row):
    for i in range(len(col)):
        parameter = col[i]
        val = row[i]
        if parameter.startswith(column_prefix["partialP"]):
            [comp, space] = parameter.split("_")
            if search_by_id(containment, space):
                add_param(
                    containment,
                    space,
                    comp[len(column_prefix["partialP"]) :],
                    float(val),
                )
            else:
                create_space(containment, space)
                add_param(
                    containment,
                    space,
                    comp[len(column_prefix["partialP"]) :],
                    float(val),
                )

        elif parameter.startswith(column_prefix["totalP"]):
            space = parameter.replace(column_prefix["totalP"], "")
            space = space[1:] if space.startswith("_") else space
            if search_by_id(containment, space):
                add_param(containment, space, "P", float(val))
            else:
                create_space(containment, space)
                add_param(containment, space, "P", float(val))


def add_temp(containment, path, time):
    file = open(path, "r")
    columns = file.readline().split()
    for line in file:
        vals = line.split()
        if vals[0] == time:
            for i in range(len(columns)):
                if columns[i].startswith("T"):
                    col = columns[i]
                    if search_by_id(containment, col[2:]):
                        add_param(containment, col[2:], "T", float(vals[i]))


def add_vol(containment, path):
    file = open(path, "r")
    for line in file:
        vals = line.split()
        if search_by_id(containment, vals[0]):
            add_param(containment, vals[0], "V", float(vals[1]))


def validate(containment):
    flg = True
    tol = 0.000001
    for space in containment:
        if space.get("T") == None:
            print("Temperature is missing in ", space.get("ID"))
            flg = False
        if space.get("P") == None:
            print("Pressure is missing in ", space.get("ID"))
            flg = False
        if space.get("V") == None:
            print("Volume is missing in ", space.get("ID"))
            flg = False
        total_x = 0
        for key in space.keys():
            if key not in not_composition:
                total_x += space.get(key)

        if total_x > 1 + tol:
            print(
                "Sum of mole fractions cannot be greater than one in ", space.get("ID")
            )
            print(total_x)
            flg = False
        elif 1 - total_x > tol:
            print(
                "Sum of mole fractions have to be one. Data missing in ",
                space.get("ID"),
            )
            flg = False
            print(1 - total_x)

    return flg


def find_data_from_file(filepath, time):
    file = open(filepath, "r")
    columns = file.readline().split()
    for line in file:
        vals = line.split()
        if vals[0] == time:
            return [columns, vals]

    return False


def process_data(CompPath, TempPath, VolPath):
    file = open(CompPath[0], "r")
    columns = file.readline().split()
    data = []

    for line in tqdm(file):
        try:
            vals = line.split()
            timestamp = vals[0]
            for i in range(len(vals)):
                vals[i] = float(vals[i])
            containment = []
            parse(containment, columns, vals)
            for filepath in CompPath:
                if filepath != CompPath[0]:
                    file_comp = find_data_from_file(filepath, timestamp)
                    if file_comp:
                        parse(containment, file_comp[0], file_comp[1])

            add_temp(containment, TempPath, timestamp)
            add_vol(containment, VolPath)

        except:
            print("Error in parsing data")

        for space in containment:
            tot = 0
            for key in space.keys():
                if key not in not_composition:
                    tot += space[key]

            for key in space.keys():
                if key not in not_composition:
                    space[key] /= tot

        if not validate(containment):
            print("Please verify the data")
            return
        res = {"TIME": vals[0]}
        for space in containment:
            if space.get("ID").upper() not in ignore_compartments:
                try:
                    ans = computeAICC(space)
                    res[space.get("ID") + "_T"] = float(ans["T"])
                    res[space.get("ID") + "_P"] = float(ans["P"])
                except:
                    print("Error in calculating AICC pressure")

        try:
            ans = homogenise(containment)
            res["Homogenised" + "_T"] = float(ans["T"])
            res["Homogenised" + "_P"] = float(ans["P"])
        except:
            print("Error in AICC pressure homogeneous mixture")
        data.append(res)

    result = pd.DataFrame.from_dict(data)
    result.to_csv("result.csv")
