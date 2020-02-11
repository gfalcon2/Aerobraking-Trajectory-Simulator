import subprocess
import time
import config
import os

def MarsGram_online(model):
    # Change input file
    keyword_month =    '  MONTH    ='
    keyword_day =      '  MDAY     ='
    keyword_year =     '  MYEAR    ='
    keyword_hour =     '  IHR      ='
    keyword_min =      '  IMIN     ='
    keyword_sec =      '  SEC      ='
    keyword_points =   '  NPOS     ='
    keyword_inlat =    '  FLAT     ='
    keyword_inlon =    '  FLON     ='
    keyword_inhgt =    '  FHGT     ='
    keyword_deltime =  '  DELTIME  ='
    keyword_dellat =   '  DELLAT   ='
    keyword_dellon =   '  DELLON   ='
    keyword_delhgt =   '  DELHGT   ='
    keyword_numberMC = '  NMONTE   ='
    keyword_seeds =    '  NR1      ='

    file_object = '/Users/giusyfalcone/MarsGram/Code/inputstd0.txt'

    lines = []
    altitude = str(model['Initial Altitude'])
    latitude = str(model['Initial Latitude'])
    longitude = str(model['Initial Longitude'])

    timereal = model['Time Real']
    year = str(timereal.year)
    month = str(timereal.month)
    day = str(timereal.day)

    hour = str(timereal.hour)
    minute = str(timereal.minute)
    sec = str(timereal.second)

    point_number = str(model['Number of Points'])
    delalt = str(model['Delta Altitude'])
    dellon = str(model['Delta Longitude'])
    dellat = str(model['Delta Latitude'])
    deltime = str(model['Delta t'])
    mc_number = str(model['Monte Carlo'])
    mc_seeds = str(int(config.index_MonteCarlo))

    with open(file_object) as infile:
        for line in infile:
            if keyword_month in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + month)
                #print(line)
            elif keyword_day in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + day)
                #print(line)
            elif keyword_year in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + year)
                #print(line)
            elif keyword_hour in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + hour)
                #print(line)
            elif keyword_min in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + minute)
                #print(line)
            elif keyword_sec in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + sec)
                #print(line)
            elif keyword_points in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + str(point_number))
                #print(line)
            elif keyword_inlat in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + latitude)
                #print(line)
            elif keyword_inlon in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + longitude)
                #print(line)
            elif keyword_inhgt in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + altitude)
                #print(line)
            elif keyword_delhgt in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + delalt)
                #print(line)
            elif keyword_dellon in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + dellon)
                #print(line)
            elif keyword_dellat in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + dellat)
                #print(line)
            elif keyword_deltime in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + deltime)
                #print(line)
            elif keyword_numberMC in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + mc_number)
                #print(line)
            elif keyword_seeds in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + mc_seeds)
            lines.append(line)

    outfile = open(file_object, 'w')
    outfile.writelines(lines)

    infile.close()
    outfile.close()

    start_time = time.time()
    args = ("/Users/giusyfalcone/MarsGram/Code/marsgram_M10.x", "-d",
            "OUTPUT.txt",)  # , "-c", "somefile.xml", "-d", "text.txt", "-r", "aString", "-f", "anotherString")
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = popen.stdout.read()
    #print("--- MARSGram execution %s seconds ---" % (time.time() - start_time))

    # Read Output File
    file_object = '/Users/giusyfalcone/MarsGram/OUTPUT.txt'
    data_list = []

    # Split Values
    with open(file_object) as fileobj:
        for line in fileobj:
            if line[5] == 'T':
                keys = line.split()
            else:
                data = line.split()
                data_list.append(data)

    # List of interesting values
    data_interesting = [[], [], [], [], [], [], [], [], [], [], []]

    # Creates lists
    for i in range(len(data_list)):
        for j in range(len(data_list[i])):
            if keys[j] == 'HgtMOLA' or keys[j] == 'HgtSFCM':
                data_interesting[0].append(float(data_list[i][j]))
            elif keys[j] == 'LatPC':
                data_interesting[1].append(float(data_list[i][j]))
            elif keys[j] == 'LonW':
                data_interesting[2].append(float(data_list[i][j]))
            elif keys[j] == 'Denkgm3':
                data_interesting[3].append(float(data_list[i][j]))
            elif keys[j] == 'Temp':
                data_interesting[4].append(float(data_list[i][j]))
            elif keys[j] == 'EWind':
                data_interesting[5].append(float(data_list[i][j]))
            elif keys[j] == 'EWTot':
                data_interesting[6].append(float(data_list[i][j]))
            elif  keys[j] == 'NWind':
                data_interesting[7].append(float(data_list[i][j]))
            elif  keys[j] == 'NWTot':
                data_interesting[8].append(float(data_list[i][j]))
            elif  keys[j] == 'VWind':
                data_interesting[9].append(float(data_list[i][j]))
            elif keys[j] == 'DensP':
                data_interesting[10].append(float(data_list[i][j]))


    list_to_remove = ['Time', 'sigD', 'Ls', 'Dust', 'LTST', 'CO2%m', 'N2%m', 'Ar%m', 'O2%m', 'CO%m', 'O%m', 'He%m',
                      'H2%m', 'H%m', 'H2O%m', 'DensP']
    for key in list_to_remove:
        keys.remove(key)


    ## MONTECARLO
    if model['Monte Carlo'] == 1:
        density = [a*b for a,b in zip(data_interesting[3],data_interesting[10])]
        data_interesting[3] = density


    config.atmospheric_data = {keys[i]: data_interesting[i] for i in range(len(keys))}

    dirpath = os.getcwd()
    test = os.listdir(dirpath)

    for item in test:
        if item.endswith(".txt"):
            os.remove(os.path.join(dirpath, item))

