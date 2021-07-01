import numpy as np
import functions as fun
from itertools import product
import copy
import csv

## Probabilidades condicionadas desde el extremo final del péptido de los aminoácidos

# Cargar matrices condicionadas 

len_aa=4

probs=[]
probarray={}
array_aa1=np.load(str("arrays/array_aa" + str(1) + ".npy"))  # Probabilidad individual del último aminoácido
array_aa1=fun.normalized_population_array(array_aa1)
probarray={'alpha' : fun.alpha_ind_prob(array_aa1), 'ppii': fun.ppii_ind_prob(array_aa1), 'beta': fun.beta_ind_prob(array_aa1),'remaining':fun.remaining_ind_prob(array_aa1)}
probs.append(probarray)
probarray={}


sec_structures=["alpha","ppii","beta","remaining"]


probs2=[]
for i in range(1,len_aa):                                                        # Probabilidades condicionadas 
    for j in sec_structures:
        probarray["alpha"]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray["ppii"]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray["beta"]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray["remaining"]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probs2.append(probarray) # 12 diccionarios, se introducen en su lista correspondiente en función del aminioácido al que se refieren
        probarray={}
    probs.append(probs2)        # 4 listas, una para cada aminoacido que contiene un diccionario dependiendo de la estructura del aminoacido vecino (0=alpha, 1=ppi, 2=beta, 3=remaining)
    probs2=[]


# Obtenemos todas las combinaciones de estructuras

sec_structures=["alpha","ppii","beta","remaining"]
combinations= (list(product(sec_structures,repeat=4)))
dictionaries=[]
dictionary={}
for i in combinations:
    dictionary={}
    dictionary[i[0]+"aa"+str(1)] = probs[0].get(str(i[0]))
    for j in range(1,len(probs)):
        if i[j-1]=="alpha":
            dictionary[i[j]+"aa"+str(j+1)]=probs[j][0].get(str(i[j]))
        if i[j-1]=="ppii":
            dictionary[i[j]+"aa"+str(j+1)]=probs[j][1].get(str(i[j]))
        if i[j-1]=="beta":
            dictionary[i[j]+"aa"+str(j+1)]=probs[j][2].get(str(i[j]))
        if i[j-1]=="remaining":
            dictionary[i[j]+"aa"+str(j+1)]=probs[j][3].get(str(i[j]))
    dictionaries.append(dictionary)

# Calculamos probabilidades de cada combinación y cambiamos la notación de las conformaciones a números y compararlos con los datos de simulación:
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

with open('probabilities/directe_conditional_probabilities.csv','w',newline='') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['structure','probabilty'])
    for row in probabilities:
        csv_out.writerow(row)     


