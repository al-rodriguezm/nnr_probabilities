import numpy as np
import functions as fun

## Importar datos

data1=np.loadtxt("rama.dat")

## Generación de matrices para cada aminoacido (matrices de 2)

data1=fun.aa_array(data1)

## Matrices de población

population_arrays=[]
for i in data1:
    population_arrays.append(fun.population_array(i))

np.save("./arrays/global_arrays.npy",population_arrays)      

# Para guardar las matrices por aminoácido.

for i in range(0,len(population_arrays)):
    name="array_aa" + str(i+1) + ".npy" 
    object= population_arrays[i]
    np.save(name,object)
