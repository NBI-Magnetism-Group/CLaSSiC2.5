import helper
import numpy as np
import os, sys


# Windows uses backward slashes for file paths and the rest uses forward slashes
if sys.platform=="win32":
    slash = "\\"
else:
    slash = "/"


def New_file(name="Saved_sim",File_number = 0):
    paths_save = slash+"CLaSSiC2.0"+slash+name+str(File_number) +".csv"
    Directory_to_here = os.getcwd()
    Complete_directory = Directory_to_here + paths_save
    if os.path.exists(Complete_directory):
        os.remove(Complete_directory)

    File = open(Complete_directory,"w")
    param, x, y, z = helper.getData()
    for key in param[0]:
        File.write(str(key)+'\n')
        File.write(str(param[0][str(key)])+'\n')

    pos = helper.getPositions()
    for i in pos:
        File.write(str(i[0])+","+str(i[1])+","+str(i[2])+"\n")

    for i in range(len(x)):
        for k in range(len(x[0,0])):
            for j in range(len(x[0])):
                File.write(str(x[i,j,k]) + ',' + str(x[i,j,k]) + ',' + str(x[i,j,k]))

                if j < len(x[0])-1:
                    File.write(',')
                
            File.write('\n')

        if i < len(x)-1:
            File.write('\n')
    File.close()
    return
    
def Append_To_File(name="Saved_sim0"):
    paths_save = slash+"CLaSSiC2.0"+slash+name+".csv"
    Directory_to_here = os.getcwd()
    Complete_directory = Directory_to_here + paths_save
    if  not os.path.exists(Complete_directory):
        print('File ' + name + ' does not exist')
        return
    
    File = open(Complete_directory,'a')
    param, x, y, z = helper.getData()
    for i in range(len(x)):
        for k in range(len(x[0,0])):
            for j in range(len(x[0])):
                File.write(str(x[i,j,k]) + ',' + str(x[i,j,k]) + ',' + str(x[i,j,k]))
                if j < len(x[0])-1:
                    File.write(',')
            File.write('\n')
        if i < len(x)-1:
            File.write('\n')

    File.close()
    return

def Load_data(name="Saved_sim0", length_of_param = 16):
    paths_save = slash+"CLaSSiC2.0"+slash+name+".csv"
    Directory_to_here = os.getcwd()
    Complete_directory = Directory_to_here + paths_save
    if  not os.path.exists(Complete_directory):
        print('File ' + name + ' does not exist')
        return
    
    File = open(Complete_directory,'r')
    data = File.read().split('\n')
    param = {}
    for i in range(length_of_param):
        param[data[i*2]] = data[i*2+1]
    
    del data[:(length_of_param)*2]
    pos = []
    for i in range(int(param['atoms'])):
        token = data[i].split(',')
        pos.append([float(token[0]),float(token[1]),float(token[2])])

    del data[:int(param['atoms'])]
    x, y, z = [], [], []
    for i in range(len(data)-1):
        spin = data[i].split(',')
        for j in range(int(len(spin)/3)):
            x.append(spin[3*j])
            y.append(spin[3*j+1])
            z.append(spin[3*j+2])

    File.close()
    return param, pos, x,y,z

def Save_last_spin(name="Saved_spin"):
    paths_save = slash+"CLaSSiC2.0"+slash+name +".dat"
    Directory_to_here = os.getcwd()
    Complete_directory = Directory_to_here + paths_save
    if os.path.exists(Complete_directory):
        os.remove(Complete_directory)

    File = open(Complete_directory,"w")
    param, x, y, z = helper.getData()
    for i in range(param[0]["atoms"]):
        sx = str(x[0,i,-1])
        sy = str(y[0,i,-1])
        sz = str(z[0,i,-1])
        if i == (param[0]["atoms"]-1):
            File.write(sx+'\n'+sy+'\n'+sz)
        else:
            File.write(sx+'\n'+sy+'\n'+sz+'\n')

    File.close()
    return
"""
This file is simply use to save the last spin orientation 
Just run it and a file will be generated
"""
Save_last_spin()
print('done')