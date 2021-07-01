import numpy as np
from numpy.lib.function_base import diff
import functions as fun
from itertools import product
import copy
import csv 


datos1=np.loadtxt("./ramaprueba.dat")
datos1=fun.aa_array(datos1)
# Función

length=4
final=[]
for a in range(0,length):
    conformations=[]
    for i in datos1[a]:
        if (i[0]>(-130) and i[0]<(-30)) and (i[1]>(-70) and i[1]<(30)):
            conformations.append("1")
        elif (i[0]>(-90) and i[0]<(-30)) and (i[1]>(115) and i[1]<(175)):
            conformations.append("2")
        elif (i[0]>(-170) and i[0]<=(-90)) and (i[1]>(120) and i[1]<(180)):
            conformations.append("3")
        else:
            conformations.append("4")
    final.append(conformations)

#

combinations=[]
for i in range(0,len(conformations)):
    combination=""
    combination=str(final[0][i]+final[1][i]+final[2][i]+final[3][i])
    combinations.append(combination)

sec_structures=["1","2","3","4"]
combinations1= (list(product(sec_structures,repeat=4)))
combinations2=[]
for i in combinations1:
    combination=""
    combination=str(i[0]+i[1]+i[2]+i[3])
    combinations2.append(combination)

probabilities=[]
for i in combinations2:
    n=0
    for j in combinations:
         if i==str(j):
            n+=1
    value=(i,n/len(datos1[0]))
    probabilities.append(value)

pos = 1
for i in range(0, len(probabilities)):
    for j in range(0, len(probabilities)-1):  
        if (probabilities[j][pos] < probabilities[j + 1][pos]):  
            temp = probabilities[j]  
            probabilities[j]= probabilities[j + 1]  
            probabilities[j + 1]= temp  

# Exportar a csv

with open('probabilities/simulation_probabilities.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['structure','probabilty'])
    for row in probabilities:
        csv_out.writerow(row)     


