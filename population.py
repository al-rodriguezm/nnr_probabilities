import numpy as np
import functions as fun
from itertools import product
import copy
import csv

# Numero total de áminoacidos del péptido y que lista queremos obtener (vecinos anteriores, posteriores o por aminoacidos (ambas))

len_aa=4

## Para poblaciones independientes de cada residuo

probarray={}
probs_ind=[]
for i in range(0,len_aa):                           # Obtenemos las probabilidades de cada una de las conformaciones en las matrices individuales, que cargamos previamente
    probarray["aa"]="independent_population"
    probarray[str("alpha_aa"+str(i+1))]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + ".npy"))))  
    probarray[str("ppii_aa"+str(i+1))]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + ".npy"))))
    probarray[str("beta_aa"+str(i+1))]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + ".npy"))))
    probarray[str("remaining_aa"+str(i+1))]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + ".npy"))))                                                    # Probabilidades condicionadas 
    probs_ind.append(probarray)
    probarray={}


## Para efecto de los residuos anteriores
# En la lista de probs: índice de 0 a len_aa-1: efecto de residuo 1 sobre 2 (hasta len_aa-1)
sec_structures=["alpha","ppii","beta","remaining"]
probarray={}
probs_pre=[]
probs2=[]
for i in range(1,len_aa):
    for j in sec_structures:
        probarray["aa"]=str("aa_"+str(i)+j)         # Obtenemos las probabilidades de cada una de las conformaciones en las matrices condicionadas, que cargamos previamente
        probarray[str("alpha_aa"+str(i+1))]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray[str("ppii_aa"+str(i+1))]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray[str("beta_aa"+str(i+1))]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probarray[str("remaining_aa"+str(i+1))]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i+1) + "_aa" + str(i) + "_" + str(j) + ".npy"))))
        probs2.append(probarray) # 12 diccionarios, se introducen en su lista correspondiente en función del aminioácido al que se refieren
        probarray={}
    probs_pre.append(probs2)        # 4 listas, una para cada aminoacido que contiene un diccionario dependiendo de la estructura del aminoacido vecino (0=alpha, 1=ppi, 2=beta, 3=remaining)
    probs2=[]


## Para efecto de los residuos posteriores
# En la lista de probs: índice de 0 a len_aa-1: efecto de residuo "len_aa" sobre "len_aa-1" (hasta len_aa-1)

sec_structures=["alpha","ppii","beta","remaining"]
probarray={}
probs_post=[]
probs2=[]
for i in reversed(range(1,len_aa)):                                                        # Probabilidades condicionadas 
    for j in sec_structures:
        probarray["aa"]=str("aa_"+str(i+1)+j)   # Obtenemos las probabilidades de cada una de las conformaciones en las matrices condicionadas, que cargamos previamente
        probarray[str("alpha_aa"+str(i))]= fun.alpha_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray[str("ppii_aa"+str(i))]= fun.ppii_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray[str("beta_aa"+str(i))]= fun.beta_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probarray[str("remaining_aa"+str(i))]= fun.remaining_ind_prob(fun.normalized_population_array(np.load(str("arrays/array_aa" + str(i) + "_aa" + str(i+1) + "_" + str(j) + ".npy"))))
        probs2.append(probarray) # 12 diccionarios, se introducen en su lista correspondiente en función del aminioácido al que se refieren
        probarray={}
    probs_post.append(probs2)        # 4 listas, una para cada aminoacido que contiene un diccionario dependiendo de la estructura del aminoacido vecino (0=alpha, 1=ppi, 2=beta, 3=remaining)
    probs2=[]


## Para generar tablas con poblacion condicionada por vecino tanto anterior como posterior

probs3=[]
probs4=[]
# Añadimos el efecto en aa1
probs4.append(probs_post[len(probs_post)-1])
for i in reversed(range(0,len(probs_post)-1)):
    probs3.append(probs_pre[len(probs_post)-2-i])
    probs3.append(probs_post[i])
    probs4.append(probs3)
    probs3=[]
# Añadimos el efecto en aa(len)
probs4.append(probs_pre[len(probs_pre)-1])


# Guardamos las tablas de población finales en un archivo csv por cada aminoácido
for i in range(0,len(probs4)):
    if i==0 or i==(len(probs4)-1):    # En el caso en el que el aminoácido solo tenga un vecino (primer y último de la cadena)
        keys=probs4[i][0].keys()
        with open(str('population/population_aa'+str(i+1)+'.csv'), 'w', newline='')  as output_file:    # Creamos un archivo con las probabilidades condicionadas
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(probs4[i])
        with open(str('population/population_aa'+str(i+1)+'.csv'), 'a', newline='') as csvfile:        # Abrimos ese nuevo archivo para introducir también las probabilidades independientes
            fieldnames = probs_ind[i].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writerow(probs_ind[i])   
    if i!=0 and i!=(len(probs4)-1):    # En el caso en el que el aminoácido tenga dos vecino
        keys=probs4[i][0][0].keys()
        with open(str('population/population_aa'+str(i+1)+'.csv'), 'w', newline='')  as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(probs4[i][0])
            dict_writer.writerows(probs4[i][1])
        with open(str('population/population_aa'+str(i+1)+'.csv'), 'a', newline='') as csvfile:
            fieldnames = probs_ind[i].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writerow(probs_ind[i])
