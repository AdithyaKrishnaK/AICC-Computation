### Notation
The state of the system is represented with a dictionary of the form
```
    r = {
        'ID':'C25',
        'P':129050,
        'T':318,
        'V':22.4,
        'H2':0.65,
        'O2':0.0735,
        'N2':0.2765
    }
```
The compositions are in mole fractions for reactants. ID is used when states for multiple compartments are required to be stored and AICC pressure and temperature calculated. A list of `[ ’P’, ’T’, ’V’, ‘ID’ ]` is stored as a variable so that when new parameters are required to be added to the state, this can be modified. All other keys except this contain the composition data. The units for pressure, temperature and volume are Pascals, Kelvin and cubic metres respectively. ID is case insensitive.

There are three main python files 

- AICC.py - To calculate AICC pressure given a composition
- conc_parser.py -  To parse data from composition data for the different compartments and find the AICC pressure and temperature for each of them
- flame_model.py - To find the rate of consumption of hydrogen, carbon monoxide and oxygen and the rate of generation of energy given a composition

