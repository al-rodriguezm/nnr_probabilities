import numpy as np
from numpy.lib.function_base import diff
import functions as fun
from itertools import product
import copy
import csv 

len_aa=4

# Cargamos todos los datos del estudio (probabilidades independientes, condicionadas empezando por ambos extremos y datos de simulación)

simulation_probs=np.loadtxt("simulacion-poblaciones-ordenados.dat",dtype=str)  # Cargamos el archivo de probabilidades calculados en la simulación
probabilities=[]

with open('probabilities/simulation_probabilities.csv', 'r') as read_obj:
    csv_reader = csv.reader(read_obj)
    simulation_probs = list(csv_reader)
    simulation_probs.pop(0)   
with open('probabilities/global_probabilities.csv', 'r') as read_obj:
    csv_reader = csv.reader(read_obj)
    global_probs = list(csv_reader)
    global_probs.pop(0)   
    probabilities.append(global_probs)
with open('probabilities/direct_conditional_probabilities.csv', 'r') as read_obj:
    csv_reader = csv.reader(read_obj)
    direct_cond_probs = list(csv_reader)  
    direct_cond_probs.pop(0)
    probabilities.append(direct_cond_probs)
with open('probabilities/inverse_conditional_probabilities.csv', 'r') as read_obj:
    csv_reader = csv.reader(read_obj)
    inverse_cond_probs = list(csv_reader)  
    inverse_cond_probs.pop(0)
    probabilities.append(inverse_cond_probs)

# Comparar todas las probabilidades con los datos de conformaciones obtenidos de la simulación

simulation_probabilities=[]
for i in simulation_probs:
    simulation_probabilities.append(list(i[0:2]))


# # Tabla comparativa de todos los valores de probabilidad

comparation=[]
for n in simulation_probabilities:
    # print(n[0])
    combination=[]
    combination.append(n[0])
    combination.append(n[1])
    
    for i in probabilities:
        j=0
        while i[j][0]!=n[0]:
            j+=1
        combination.append(i[j][1])
    comparation.append(combination)

# Exportar csv

with open('probabilities/combination_probabilities.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['combination','simulation','global_prob','direct_cond_prob','inverse_cond_prob'])
    for row in comparation:
        csv_out.writerow(row)     

# Diferencias
# Para calcular las diferencias entre los datos de simulaciones y las aproximaciones realizadas (independent y NNR)

differences=[]
for n in simulation_probabilities:
    # print(n[0])
    combination=[]
    combination.append(n[0])
    combination.append(n[1])

    for i in probabilities:
        j=0
        while i[j][0]!=n[0]:
            j+=1
        combination.append(abs(float(i[j][1])-float(n[1])))
    differences.append(combination)

# Exportar csv

with open('probabilities/combination_diferences.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['combination','simulation','global_prob','direct_cond_prob','inverse_cond_prob'])
    for row in differences:
        csv_out.writerow(row)     
