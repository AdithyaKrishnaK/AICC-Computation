from math import exp, pow

not_composition = ["P", "V", "T", "ID"]


def flame(r, data={"F": 0.541, "x_c": 0.07, "x_o2": 0.05, "x_d": 0.55}):
    R = 8.314
    n_r = r["P"] * r["V"] / (R * r["T"])

    if r.get("H2") == None:
        r["H2"] = 0

    if r.get("O2") == None:
        r["O2"] = 0

    if r.get("H2O") == None:
        r["H2O"] = 0

    if r.get("CO") == None:
        r["CO"] = 0

    # Calculate Burning Velocity
    x = r["H2"] + r["CO"]
    if x > 0.42:
        a1 = 4.644e-4
        a2 = 9.898e-4
        a3 = -1.264e-3
        a4 = 1.571
        a5 = -2.476e-1
        a6 = -2.24
    else:
        a1 = 4.644e-4
        a2 = -2.119e-3
        a3 = 2.344e-3
        a4 = 1.571
        a5 = 3.839e-1
        a6 = -2.21
    B = a1 + a2 * (0.42 - x) + a3 * (0.42 - x) * (0.42 - x)
    C = a4 + a5 * (0.42 - x)
    D = a6
    vf = B * pow(r["T"], C) * exp(D * r["H2O"])
    l = pow(r["V"], 0.3333)
    deltaT = l / vf
    rate = {}

    # Parameters to check if explosion happens
    x = r["H2"] + data["F"] * r["CO"]

    x_d = 0
    for key in r.keys():
        if key in ["CO2", "D2O", "H2O"]:
            x_d += r[key]

    if r.get("N2") != None and r.get("O2") != None:
        x_d += 0.79 * (r["N2"] - 3.774 * r["O2"])

    # checking if explosion happens
    if x >= data["x_c"] and r["O2"] >= data["x_o2"] and x_d <= data["x_d"]:
        if (r["H2"] + r["CO"]) / 2 > r["O2"]:
            # Less oxygen
            rate["O2"] = -r["O2"] * n_r / deltaT
            rate["H2"] = rate["O2"]
            rate["CO"] = rate["O2"]
        else:
            # Excess oxygen
            rate["H2"] = -r["H2"] * n_r / deltaT
            rate["CO"] = -r["CO"] * n_r / deltaT
            rate["O2"] = 0.5 * (rate["H2"] + rate["CO"])
    else:
        print("No explosion")
        return False

    rate["H2O"] = -rate["H2"]
    rate["CO2"] = -rate["CO"]
    # Enthalpy of combustion
    q_h2 = 285.8  # KJ/mol
    q_co = 283.0  # kJ/mol
    energy = -(q_h2 * rate["H2"] + q_co * rate["CO"])
    return [rate, energy]
