def missiondef(mission):
    words = ""
    sentence = mission['Purpose']
    wordsinpurpose = []

    for letter in sentence:

        if (letter == ' ') or (letter == '.'):
            wordsinpurpose.append(words)
            words = ""
        else:
            words += str(letter)

    [e, d, l, a] = [0, 0, 0, 0]

    for words in wordsinpurpose:
        if (words == 'entry') or (words == 'Entry'):
            e = 1
            #print(words)
        if (words == 'descent') or (words == 'Descent'):
            d = 1
            #print(words)
        if (words == 'landing') or (words == 'Landing'):
            l = 1
            #print(words)
        if (words == 'aerobraking') or (words == 'Aerobraking'):
            a = 1
            #print(words)
        if (words == 'Earth'):
            p = 0
            #print(words)
        elif (words == 'Mars'):
            p = 1
            #print(words)
        elif (words == 'Venus'):
            p = 2
            #print(words)

    class Mission:
        def __init__(self, e, d, l, a, p):
            self.e = e
            self.d = d
            self.l = l
            self.a = a
            self.planet = p

    M = Mission(e,d,l,a,p)

#Gravity Model
    if (mission['Gravity Model'] == 'constant') or (mission['Gravity Model'] == 'Constant'):
        gm = 0
        #print('Gravity Model = Constant')
    elif (mission['Gravity Model'] == 'Inverse squared') or (mission['Gravity Model'] == 'inverse squared') or (mission['Gravity Model'] == 'Inverse Squared') :
        gm = 1
        #print('Gravity Model = Inverse squared')
    elif (mission['Gravity Model'] == 'Inverse Squared and J2 effect') or (mission['Gravity Model'] == 'inverse squared and J2 effect') or (mission['Gravity Model'] == 'Inverse quared and J2 effect'):
        gm = 2
        #print('Gravity Model = Inverse Squared and J2 effect')
    else:
        gm = 1
        #print('Gravity Model = Inverse squared')



#Density Model
    if (mission['Density Model'] == 'constant') or (mission['Density Model'] == 'Constant'):
        dm = 0
        #print('Density Model = Constant')
    elif (mission['Density Model'] == 'exponential') or (mission['Density Model'] == 'Exponential'):
        dm = 1
        #print('Density Model = Exponential')
    elif (mission['Density Model'] == 'No-Density') or (mission['Density Model'] == 'No-density'):
        dm = 2
        #print('Density Model = No-Density')
    elif (mission['Density Model'] == 'MARSGram') or (mission['Density Model'] == 'MarsGram'):
        dm = 3
        #print('Density Model = Mars Gram')
    else:
        dm = 1
        #print('Density Model = Exponential')

# MonteCarlo
    wm = int(mission['Wind'])

# Aerodynamic Model
    if (mission['Aerodynamic Model'] == 'Cd and Cl Constant') or (mission['Aerodynamic Model'] == 'Cd and Cl constant'):
        am = 0
        #print('Aerodynamic Model = Cd and Cl Constant')
    elif (mission['Aerodynamic Model'] == 'Diffusive') or (mission['Aerodynamic Model'] == 'Mach-dependent'):
        am = 1
        #print('Aerodynamic Model = Mach-dependent')
    elif (mission['Aerodynamic Model'] == 'No-Ballistic flight with axial coefficient') or (mission['Aerodynamic Model'] == 'No-ballistic flight with axial coefficient'):
        am = 2
        #print('Aerodynamic Model = No-Ballistic flight')
    else:
        am = 0
        #print('Aerodynamic Model = Ballistic flight')

# Control Model
    if (mission['Control'] == 1):
        cm = 1
    elif (mission['Control'] == 0):
        cm = 0
    else:
        cm = 0

# Thermal Model
    words = ""
    sentence = mission['Thermal Model']
    wordsinpurpose = []

    for letter in sentence:
        if letter != ' ':
            words += str(letter)
        else:
            wordsinpurpose.append(words)
            words = ""

    [c, r, t] = [0, 0, 0]

    for words in wordsinpurpose:
        if (words == 'convective') or (words == 'Convective'):
            c = 1
        if (words == 'radiative') or (words == 'Radiative'):
            r = 1
        if (words == 'Maxwellian') or (words == 'Shaaf and Chambre'):
            t = 1

    class ThermalModel:
        def __init__(self, c, r, t):
            self.c = c
            self.r = r
            self.t = t

    tm = ThermalModel(c,r,t)
    #if (tm.c == 1 ) and (tm.r == 1):
        #print('Thermal Model = Convective Heating (Sutton-Graves) + Radiative Heating')
    #elif (tm.c == 1) and (tm.r == 0):
        #print('Thermal Model = Convective Heating (Sutton-Graves)')
    #elif (tm.c == 0) and (tm.r == 1):
        #print('Thermal Model = Radiative Heating (Sutton-Graves)')
    #elif (tm.c == 0) and (tm.r == 0) and (tm.t == 1):
        #print('Thermal Model = Maxwellian Model for Flat Plate')

# MonteCarlo
    mc = int(mission['Monte Carlo'])

    class InitialParameters:
        def __init__(self, M, gm, dm, wm, am, tm, cm, mc):
            self.M= M
            self.gm = gm
            self.dm = dm
            self.wm = wm
            self.am = am
            self.tm = tm
            self.cm = cm
            self.mc = mc
    global ip
    ip = InitialParameters(M,gm,dm,wm,am,tm,cm,mc)
    return ip
