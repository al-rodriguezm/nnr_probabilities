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

with open('probabilities/independent_probabilities.csv', 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj)
    # Pass reader object to list() to get a list of lists
    independent_probs = list(csv_reader)
    independent_probs.pop(0)   
    probabilities.append(independent_probs)
with open('probabilities/directe_conditional_probabilities.csv', 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj)
    # Pass reader object to list() to get a list of lists
    directe_cond_probs = list(csv_reader)  
    directe_cond_probs.pop(0)
    probabilities.append(directe_cond_probs)
with open('probabilities/inverse_conditional_probabilities.csv', 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj)
    # Pass reader object to list() to get a list of lists
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

with open('probabilities/probabilities_comparation.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['combination','indepent_prob','directe_cond_prob','inverse_cond_prob'])
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

for i in differences:
    print(i)


# Exportar csv

with open('probabilities/differences_comparation.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['combination','indepent_prob','directe_cond_prob','inverse_cond_prob'])
    for row in differences:
        csv_out.writerow(row)     
