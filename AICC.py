from functools import partial
from scipy import optimize
import pyromat as pm


not_composition = ["P", "V", "T", "ID"]


def uchange(r, p, n_p, t):
    pm.config["unit_matter"] = "kmol"
    data = {}
    R = 8.314
    n_r = r["P"] * r["V"] / (R * r["T"])
    for key in r.keys():
        if key not in not_composition:
            data[key] = pm.get("ig." + key)
    for key in p.keys():
        if key not in not_composition and data.get(key) == None:
            data[key] = pm.get("ig." + key)
    val = 0
    for key in r.keys():
        if key not in not_composition:
            val -= data[key].e(T=r["T"]) * r[key] * n_r

    for key in p.keys():
        if key not in not_composition:
            val += data[key].e(T=t) * p[key] * n_p
    return val

#data contains
# F - Scaling factor to find equivalent hydrogen for carbon monoxide
# x_c - Minimum mole fraction of equivalent hydrogen
# x_o2 - Minimum mole fraction of oxygen
# x_d - Maximum mole fraction of dilutents

def computeProd(r, data):
    p = {}
    R = 8.314
    n_r = r["P"] * r["V"] / (R * r["T"])

    #Checking for hydrogen or deutrium
    if r.get("H2") != None:
        h_key = "H2"
    else:
        h_key = "D2"

    # Setting values zero so that later on error is not raised
    if r.get("CO") == None:
        r["CO"] = 0

    if r.get(h_key + "O") == None:
        r[h_key + "O"] = 0

    if r.get("CO2") == None:
        r["CO2"] = 0

    x = r[h_key] + data["F"] * r["CO"]

    x_d = 0
    for key in r.keys():
        if key in ["CO2", "D2O", "H2O"]:
            x_d += r[key]

    if r.get("N2") != None and r.get("O2") != None:
        x_d += 0.79 * (r["N2"] - 3.774 * r["O2"])

    if x >= data["x_c"] and r["O2"] >= data["x_o2"] and x_d <= data["x_d"]:
        if (r[h_key] + r["CO"]) / 2 > r["O2"]:
            # Oxygen limted - distribution according to proportion
            k = 2 * r["O2"] / (r[h_key] + r["CO"])
            p[h_key + "O"] = r[h_key + "O"] * n_r + k * r[h_key] * n_r
            p["CO2"] = r["CO2"] * n_r + k * r["CO"] * n_r
            p["CO"] = (r["CO"] - k * r["CO"]) * n_r
            p[h_key] = (r[h_key] - k * r[h_key]) * n_r
            if p[h_key] < 0:
                print("Error: Concentration cannot be negative")
                p[h_key] = 0
            if p["CO"] < 0:
                print("Error: Concentration cannot be negative")
                p["CO"] = 0
            p["O2"] = 0
        else:
            #excess oxygen
            p[h_key + "O"] = r[h_key + "O"] * n_r + r[h_key] * n_r
            p[h_key] = 0
            p["CO"] = 0
            p["CO2"] = r["CO2"] * n_r + r["CO"] * n_r
            p["O2"] = (r["O2"] - (r[h_key] + r["CO"]) / 2) * n_r
    else:
        # No explosion
        return r

    # No changes to non reacting species
    for key in r.keys():
        if key not in list(
            set().union(["H2", "D2", "O2", "D2O", "H2O", "CO", "CO2"], not_composition)
        ):
            p[key] = n_r * r[key]
    return p


def inbuiltsolver(fun, xg, x1):
    try:
        sol = optimize.newton(fun, xg, x1=x1)
    except:
        print("No convergence")
        return 0

    if abs(fun(sol)) > 0.001 or sol < 273:
        print("Error: Wrong convergence")
        return 0
    return sol

#data contains
# F - Scaling factor to find equivalent hydrogen for carbon monoxide
# x_c - Minimum mole fraction of equivalent hydrogen
# x_o2 - Minimum mole fraction of oxygen
# x_d - Maximum mole fraction of dilutents
def computeAICC(r, data={"F": 0.541, "x_c": 0.07, "x_o2": 0.05, "x_d": 0.55}):

    if r.get("H2") == None and r.get("D2") == None:
        r["H2"] = 0

    if r.get("H2") != None and r.get("D2") != None:
        print("Error, both hydrogen and deutrium present. Not Supported")
        return r

    if r.get("O2") == None:
        print("No oxygen. No explosion")
        return r

    R = 8.314
    p = computeProd(r, data)

    # No reaction
    if p.get("T") != None:
        return p

    n_p = 0
    for key in p.keys():
        if key not in not_composition:
            n_p += p[key]

    # Converting product composition to mole fraction
    for key in p.keys():
        if key not in not_composition:
            p[key] /= n_p

    # Second guess for Newton-Raphson method
    xg = float(-uchange(r, p, n_p, r["T"]) / (2.5 * R * n_p) + r["T"])
    if xg > 6000:
        xg = 6000

    fun = partial(uchange, r, p, n_p)
    T = inbuiltsolver(fun, r["T"], xg)
    P_AICC = n_p * R * T / r["V"]
    p["T"] = T
    p["P"] = P_AICC
    return p
