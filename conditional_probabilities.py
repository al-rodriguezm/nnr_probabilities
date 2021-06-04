import numpy as np
import functions as fun
from itertools import product
import copy

# Cargar todas las matrices condicionadas.

# Inversas

len_aa=4

probs=[]
probarray={}
array_aa1=np.load(str("array_aa" + str(len_aa-1) + ".npy"))  # Probabilidad individual del último aminoácido
array_aa1=fun.normalized_population_array(array_aa1)
probarray={'alpha' : fun.alpha_ind_prob(array_aa1), 'ppii': fun.ppii_ind_prob(array_aa1), 'beta': fun.beta_ind_prob(array_aa1),'remaining':fun.remaining_ind_prob(array_aa1)}
probs.append(probarray)
probarray={}

# sec_structures=["alpha","ppii","beta","remaining"]

# for i in range(1,len_aa):                                                        # Probabilidades condicionadas 
#     for j in sec_structures:
#         probarray[str("array_aa"+str(i)+ "_alpha" + "_aa"+str(i+1)+"_"+str(j))]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
#         probarray[str("array_aa"+str(i)+ "_ppii" + "_aa"+str(i+1)+"_"+str(j))]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
#         probarray[str("array_aa"+str(i)+ "_beta" + "_aa"+str(i+1)+"_"+str(j))]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
#         probarray[str("array_aa"+str(i)+ "_remaining" + "_aa"+str(i+1)+"_"+str(j))]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
#         # probs.append(probarray) 12 diccionarios, uno para cada estructura de cada aminoacido
#         # probarray={}
#     probs.append(probarray)        # 4 diccionarios, uno para cada aminoácido
#     probarray={}
# print(probs[1])

sec_structures=["alpha","ppii","beta","remaining"]


probs2=[]
for i in range(1,len_aa):                                                        # Probabilidades condicionadas 
    for j in sec_structures:
        probarray["alpha"]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray["ppii"]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray["beta"]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray["remaining"]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probs2.append(probarray) # 12 diccionarios, se introducen en su lista correspondiente en función del aminioácido al que se refieren
        probarray={}
    probs.append(probs2)        # 4 listas, una para cada aminoacido que contiene un diccionario dependiendo de la estructura del aminoacido vecino (0=alpha, 1=ppi, 2=beta, 3=remaining)
    probs2=[]
    # probarray={}

# print(probs)

# print("######")

sec_structures=["alpha","ppii","beta","remaining"]
combinations= (list(product(sec_structures,repeat=4)))
# dictionaries=[]
# for i in combinations:
#     for j in range(1,len(i)):
#         dictionary={} 
#         dictionary[i[j]+"aa"+str(1)] =  probs[0][j].get(str(i[j]))
#         dictionaries.append(dictionary)
#     if i[j] == "alpha":
#         dictionary={} 
#         dictionary[i[j]+"aa"+str(2)] =  probs[1][0].get(str(i[j]))
#         dictionaries.append(dictionary)
#     if i[j] == "ppii":
#         dictionary={} 
#         dictionary[i[j]+"aa"+str(2)] =  probs[1][1].get(str(i[j]))
#         dictionaries.append(dictionary)
#     if i[j] == "beta":
#         dictionary={} 
#         dictionary[i[j]+"aa"+str(2)] =  probs[1][2].get(str(i[j]))
#         dictionaries.append(dictionary)
#     if i[j] == "remaining":
#         dictionary={} 
#         dictionary[i[j]+"aa"+str(2)] =  probs[1][3].get(str(i[j]))
#         dictionaries.append(dictionary)


dictionaries=[]
dictionary={}
for i in combinations:
    for j in range(1,len(probs)):
        if i[j-1]=="alpha":
            dictionary[i[j]+"aa"+str(j)]=probs[j][0].get(str(i[j]))
        if i[j-1]=="ppii":
            dictionary[i[j]+"aa"+str(j)]=probs[j][1].get(str(i[j]))
        if i[j-1]=="beta":
            dictionary[i[j]+"aa"+str(j)]=probs[j][2].get(str(i[j]))
        if i[j-1]=="remaining":
            dictionary[i[j]+"aa"+str(j)]=probs[j][3].get(str(i[j]))
    dictionary[i[0]+"aa"+str(4)] = probs[0].get(str(i[0]))
    dictionary={}
    dictionaries.append(dictionary)


# Directas
# dictionaries=[]
# dictionary={}
# for i in combinations:
#     dictionary={}
#     dictionary[i[0]+"aa"+str(1)] = probs[0].get(str(i[0]))
#     for j in range(1,len(probs)):
#         if i[j-1]=="alpha":
#             dictionary[i[j]+"aa"+str(j+1)]=probs[j][0].get(str(i[j]))
#         if i[j-1]=="ppii":
#             dictionary[i[j]+"aa"+str(j+1)]=probs[j][1].get(str(i[j]))
#         if i[j-1]=="beta":
#             dictionary[i[j]+"aa"+str(j+1)]=probs[j][2].get(str(i[j]))
#         if i[j-1]=="remaining":
#             dictionary[i[j]+"aa"+str(j+1)]=probs[j][3].get(str(i[j]))
#     dictionaries.append(dictionary)
# # n=0
# for i in dictionaries:
#       print(i)

# print(n)

probabilities=[]
for i in dictionaries:
    values=list(i.values())
    keys=list(i.keys())
    total_probability=1
    peptide=""
    value=()
    for j in range(0,len(i)):
        total_probability*=values[j]
        if j==0:
            peptide+=keys[j]
        else: 
            peptide+="-"+keys[j]
    
    value=(peptide,total_probability)
    probabilities.append(value)
    
    values=[]
    keys=[]


# print(probabilities)

# # Ordenamos el diccionario

pos = 1
for i in range(0, len(probabilities)):
    for j in range(0, len(probabilities)-1):  
        if (probabilities[j][pos] < probabilities[j + 1][pos]):  
            temp = probabilities[j]  
            probabilities[j]= probabilities[j + 1]  
            probabilities[j + 1]= temp  

# Las ocho combinaciones de estructuras más probables

for i in probabilities[0:8]:
    print(i)

# print(len(dictionaries))

# print(dictionaries)
# probabilities=[]
# for i in dictionaries:
#     values=list(i.values())
#     keys=list(i.keys())
#     total_probability=1
#     peptide=""
#     value=()
#     for j in range(0,len(i)):
#         total_probability*=values[j]
#         if j==0:
#             peptide+=keys[j]
#         else: 
#             peptide+="-"+keys[j]
    
#     value=(peptide,total_probability)
#     probabilities.append(value)
    
#     values=[]
#     keys=[]