import numpy as np
import functions as fun
import itertools
import copy

## Importar datos

datos1=np.loadtxt("prueba.dat")

## Matrices de población (total)


max=len((datos1[1])/2)-1

for i in range(0,int(max)):
    if i==0:
        name = "arrays/array_" + "a" + str(i) + "_" + str(i+1) 
        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i+1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i+1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i+1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i+1))
        name = ""

    if i==max:
        name = "arrays/array_" + "a" + str(max)= + "_" + str(max+1) 
        
        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i-1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i-1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i-1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i-1))
        name = ""

    else:
        name = "arrays/array_" + "a" + str(i)= + "_" + str(i+1)

        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i+1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i+1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i+1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i+1))
        name = ""
        
        name = "arrays/array_" + "a" + str(i-1)= + "_" + str(i-1) 

        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i-1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i-1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i-1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i-1))
        name = ""


## Condicional a partir de aminoácido inicial

max=len((datos1[1])/2)-1

for i in range(0,int(max)):

    if i==max:
        name = "arrays/array_" + "a" + str(max)= + "_" + str(max+1) 
        
        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i-1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i-1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i-1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i-1))
        name = ""

    else:
        name = "arrays/array_" + "a" + str(i)= + "_" + str(i+1)

        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i+1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i+1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i+1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i+1))
        name = ""



## Condicional a partir de aminóacido final

max=len((datos1[1])/2)-1

for i in range(0,int(max)):
    if i==0:
        name = "arrays/array_" + "a" + str(i) + "_" + str(i+1) 
        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i+1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i+1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i+1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i+1))
        name = ""

    else:
        
        name = "arrays/array_" + "a" + str(i-1)= + "_" + str(i-1) 

        np.save(str(name + "alpha"),fun.population_conditional_array_alpha(datos1, i, i-1))
        np.save(str(name + "ppii"), fun.population_conditional_array_ppii(datos1, i, i-1))
        np.save(str(name + "beta"), fun.population_conditional_array_beta(datos1, i, i-1))
        np.save(str(name + "remaining"),  fun.population_conditional_array_remaining(datos1, i, i-1))
        name = ""
