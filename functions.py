import numpy as np
import copy 

## Generaricón de matrices para cada aminoacido (matrices de 2)

def aa_array(array):
    """
    (numpy.ndarray) -> list
    Extrae de los datos de rama.dat las matrices con los valores de angulos diedros a lo largo de la simulación para cada uno de los aminoácidos que conforman el péptido.    
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de los aminóacidos del péptido en cuestión.
        Salida: Una lista conformada con las matrices de datos de angulos diedros de cada aminoácido en la simulación. 
    """
    k=0                     # Obtenemos el número de columnas del array para conocer el número de aminoácidos del péptido
    while k==0:
        for i in array:
            k=len(i)

    i=0                    # Dividimos la array general en matrices de dos columnas con los datos de ángulos phi y psi respectivamente y las introducimos en una lista
    j=2
    arrays=[]
    while j<=k:
        aa_array=array[:,i:j]
        arrays.append(aa_array)
        i=i+2
        j=j+2

    return(arrays)

## Transformacion a array de 0 a 360

def population_array(array):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Obtiene de una array con los valores de ángulos diedros para un aminoácido, una array de un tamaño 360x360 con la población de cada uno de los puntos, los cuales 
    representan una combinación de ángulos diedros específica. 
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de un aminóacidos de la cadena.
        Salida: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente.
    """

    array[:,0]=array[:,0]+180
    array[:,1]=-1*array[:,1]+180
    print(array)
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi
        k=0                                     # Indica que se han encontrado los ángulos correspondientes
        while (a[0]>i or a[1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
            # print(a[0],i)
            # print(a[1],j)
            if (a[0]<=i and a[1]<=j):           # i>360 -> matriz_vacia[0,j]        j>360 -> matriz_vacia[i,0] o bien i-1 j -1 (para 0-359)
                empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[0]>i:
                    i=i+1
                if a[1]>j:
                    j=j+1
    return empty_array


## Normalizacion de matrices

def normalized_population_array(array):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Normaliza la array de población dividiendo cada valor entre el número total de puntos, de manera que la suma total de los valores de la array sea 1.
        Entrada: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente
        Salida: La misma array con valores normalizados a 1.
    """
    
    n=0                 # Obtenemos el valor de la poblacion total (numero de filas)
    s=0
    for i in array: 
        for j in i:
            if j!=0:
                n=n+j
                
    population_array=[]       # Dividimos los valores de la array entre la poblacion total para normalizar los datos
    for i in array:
        population_array.append(i/n)
    return(population_array)

def check(array):       # Para comprobar que se han tenido en cuenta todos los valores de la array
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Normaliza la array de población dividiendo cada valor entre el número total de puntos, de manera que la suma total de los valores de la array sea 1.
        Entrada: array con valores normalizados
        Salida: Número total de puntos en la array
    """
    n=0
    for i in array:
        for j in i:
            if j!=0:
                n=n+j
    return(n)

## Probabilidades por NNR
def population_conditional_array_alpha(array,aa1,aa2):
    array1=copy.deepcopy(array)
    array1[:,2*aa1]=array1[:,2*aa1]+180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los �ngulos correspondientes
        k=0
        # print(a[2*aa2])
        # print(a[2*aa2+1])
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de �ngulos psi y phi de cada fila
            # print(a[0],i)
            # print(a[1],j)
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):   
                if (a[2*aa2]>= -130 and a[2*aa2]<= -30) and (a[2*aa2+1]>= -70 and a[2*aa2+1]<=30):  # Incluimos la restricci�n de conformaci�n alfa en el amino�cido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1

    return empty_array

# print(sum(itertools.chain.from_iterable(population_conditional_array_alpha(datos1, 1, 0))))  # Matriz de aa 1 condicionado por aa 2 sea beta
# print(population_conditional_array_alpha(datos1, 1, 0)[98][41])


def population_conditional_array_ppii(array,aa1,aa2):
    array1=copy.deepcopy(array)
    array1[:,2*aa1]=array1[:,2*aa1]+180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los �ngulos correspondientes
        k=0
        # print(a[2*aa2])
        # print(a[2*aa2+1])
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de �ngulos psi y phi de cada fila
            # print(a[0],i)
            # print(a[1],j)
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):   
                if (a[2*aa2]>= -90 and a[2*aa2]<= -30) and (a[2*aa2+1]>= 115 and a[2*aa2+1]<=175):  # Incluimos la restricci�n de conformaci�n ppii en el amino�cido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1
    return empty_array


# print(sum(itertools.chain.from_iterable(population_conditional_array_ppii(datos1, 0, 1))))  # Matriz de aa 1 condicionado por aa 2 sea beta
# print(population_conditional_array_ppii(datos1, 0, 1)[99][200])


def population_conditional_array_beta(array,aa1,aa2):
    array1=copy.deepcopy(array)
    array1[:,2*aa1]=array1[:,2*aa1]+180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los �ngulos correspondientes
        k=0
        # print(a[2*aa2])
        # print(a[2*aa2+1])
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de �ngulos psi y phi de cada fila
            # print(a[0],i)
            # print(a[1],j)
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):   
                if (a[2*aa2]>= -170 and a[2*aa2]<= -90) and (a[2*aa2+1]>= 120 and a[2*aa2+1]<=180):  # Incluimos la restricci�n de conformaci�n beta en el amino�cido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1
    return empty_array

# print(sum(itertools.chain.from_iterable(population_conditional_array_beta(datos1, 1, 2))))  # Matriz de aa 1 condicionado por aa 2 sea beta
# print(population_conditional_array_beta(datos1, 2, 3)[99][200])

def population_conditional_array_remaining(array,aa1,aa2):
    array1=copy.deepcopy(array)
    array1[:,2*aa1]=array1[:,2*aa1]+180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los �ngulos correspondientes
        k=0
        # print(a[2*aa2])
        # print(a[2*aa2+1])
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de �ngulos psi y phi de cada fila
            # print(a[0],i)
            # print(a[1],j)
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):   
                if not(a[2*aa2]>= -130 and a[2*aa2]<= -30) or not(a[2*aa2+1]>= -70 and a[2*aa2+1]<=30):       # No alpha
                    if not(a[2*aa2]>= -90 and a[2*aa2]<= -30) or not(a[2*aa2+1]>= 115 and a[2*aa2+1]<=175):  # No ppii
                        if not(a[2*aa2]>= -170 and a[2*aa2]<= -90) or not(a[2*aa2+1]>= 120 and a[2*aa2+1]<=180):  # No beta
                            empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1

    array[:,2*aa1]=array[:,2*aa1]-180
    array[:,2*aa1+1]=1*array[:,2*aa1+1]+180
    return empty_array

def alpha_ind_prob(array):
    phi=-80+180
    psi=20+180
    n=0
    for i in range(phi-50,phi+50):
        for j in range(psi-50,psi+50):
            n=n+array[i-1][j-1]
    return(n)

def ppii_ind_prob(array):
    phi=-60+180
    psi=-145+180
    n=0
    for i in range(phi-30,phi+30):
        for j in range(psi-30,psi+30):
            n=n+array[i-1][j-1]
    return(n)

def beta_ind_prob(array):
    phi=-130+180
    psi=-140+180
    n=0
    for i in range(phi-40,phi+40):
        for j in range(psi-40,psi+40):
            n=n+array[i-1][j-1]
    return(n)

def remaining_ind_prob(array):
    n=1-(alpha_ind_prob(array)+ppii_ind_prob(array)+beta_ind_prob(array))
    return(n)
    
# print(sum(itertools.chain.from_iterable(population_conditional_array_remaining(datos1, 0, 1))))  # Matriz de aa 1 condicionado por aa 2 sea beta
# print(population_conditional_array_remaining(datos1, 0, 1)[99][200])

def main():
    print("Módulo con las funciones utilizadas para el transcurso del proyecto de Trabajo de Fin de Máster")

if __name__ == "__main__":
    main()