import numpy as np
import functions as fun

## Importar datos

datos1=np.loadtxt("rama.dat")

## Generaricón de matrices para cada aminoacido (matrices de 2)

datos1=fun.aa_array(datos1)

## Matrices de población

population_arrays=[]
for i in datos1:
    population_arrays.append(fun.population_array(i))

np.save("./arrays/independent_arrays.npy",population_arrays)      

# Para guardar las matrices por aminoácido.

for i in range(0,len(population_arrays)):
    name="array_aa" + str(i+1) + ".npy" 
    object= population_arrays[i]
    np.save(name,object)

## Normalizacion de matrices

normalized_population_arrays=[]
for i in population_arrays:
    normalized_population_arrays.append(fun.normalized_population_array(i))