import numpy as np
import functions as fun

## Conversion de las matrices .npy de interés en .txt para representarlas en los scripts R

path=str(input("Select an array:"))

array=np.load(path)

normalized_array=fun.normalized_population_array(array)

if path[len(path)-3:]=="npy":
    new_name=path[:-3]+"txt"
    np.savetxt(new_name,normalized_array)

print(new_name)