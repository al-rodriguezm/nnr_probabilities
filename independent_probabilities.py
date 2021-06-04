import numpy as np
from itertools import product
import functions as fun
import csv

## Probabilidades independientes de los aminoácidos

len_aa=4 

array=np.load("arrays/independent_arrays.npy")

normalized_population_arrays=[]
for i in array:
    normalized_population_arrays.append(fun.normalized_population_array(i))

# Obtenemos las probabilidades de la estructura para cada aminoacido del peptido 
probs=[]
probarray={}
for i in normalized_population_arrays:
    probarray={'alpha' : fun.alpha_ind_prob(i), 'ppii': fun.ppii_ind_prob(i), 'beta': fun.beta_ind_prob(i),'remaining':fun.remaining_ind_prob(i)}
    probs.append(probarray)
    probarray={}

for i in probs:
     print(i)

sec_structures=["alpha","ppii","beta","remaining"]
combinations= (list(product(sec_structures,repeat=4)))
dictionaries=[]
for i in combinations:
    dictionary={}
    for j in range(0,len(probs)): 
        dictionary[i[j]+"aa"+str(j+1)] =  probs[j].get(str(i[j]))
    dictionaries.append(dictionary)

# Calculamos probabilidades de cada combinación yy cambiamos la notación de las conformaciones a números y compararlos con los datos de simulación:
# 1. alpha
# 2. ppII
# 3. beta
# 4. remaining

probabilities=[]
for i in dictionaries:
    values=list(i.values())
    keys=list(i.keys())
    total_probability=1
    peptide=""
    value=()
    for j in range(0,len(i)):
        total_probability*=values[j]
        if keys[j]==str("alphaaa"+str(j+1)):
            peptide+="1"
        if keys[j]==str("ppiiaa"+str(j+1)):
            peptide+="2"
        if keys[j]==str("betaaa"+str(j+1)):
            peptide+="3"
        if keys[j]==str("remainingaa"+str(j+1)):
            peptide+="4"
    
    value=(peptide,total_probability)
    probabilities.append(value)
    
    values=[]
    keys=[]
    
# Ordenamos el diccionario

pos = 1
for i in range(0, len(probabilities)):
    for j in range(0, len(probabilities)-1):  
        if (probabilities[j][pos] < probabilities[j + 1][pos]):  
            temp = probabilities[j]  
            probabilities[j]= probabilities[j + 1]  
            probabilities[j + 1]= temp  

# Exportar a csv

with open('probabilities/independent_probabilities.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['structure','probabilty'])
    for row in probabilities:
        csv_out.writerow(row)     



