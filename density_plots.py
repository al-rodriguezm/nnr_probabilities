import numpy as np
import functions as fun
from itertools import product
import copy
import csv

# Elegimos la matriz que queremos representar 

import numpy as np
import functions as fun

path=str(input("Select an array:"))

array=np.load(path)

normalized_array=fun.normalized_population_array(array)

if path[len(path)-3:]=="npy":
    new_name=path[:-3]+"txt"
    np.savetxt(new_name,normalized_array)



# # Para obtener poblaciones por phi y psi (Util para herramienta Ramachandran si la usamos)

# densities=[]
# densities.append(["phi","psi","density"])
# row=[]
# for i in range(0,360):
#     for j in range(0,360):
#         row.append(i-180)
#         row.append(-j-180)
#         row.append(array[i][j])
#         densities.append(row)
#         row=[]

# print(densities[0:14])

# with open("plots/array_aa2_aa1_beta_plot.csv", "w", newline="") as output_file:
#     writer = csv.writer(output_file)
#     writer.writerows(densities)

# # Para diagram de Ramachandran de todo el péptido, sumamos todas las poblaciones de las matrices de aminoácidos independientes.

# array=np.load("arrays/independent_arrays.npy")
# empty_array=np.zeros((360,360),dtype=int)
# for i in array:
#     empty_array=empty_array+fun.normalized_population_array(i)

# np.savetxt("plots/ramachandran_plot.txt",empty_array)