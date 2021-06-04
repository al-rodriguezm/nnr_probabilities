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
    array[:,0]=array[:,0]+180                  # La primera fila de la matriz representa el ángulo -180
    array[:,1]=-1*array[:,1]+180               # La primera columna de la matriz representa el ángulo 180
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array:
        i=0                                     # Ángulo phi
        j=0                                     # Ángulo psi
        k=0                                     # Indica que se han encontrado los ángulos correspondientes
        while (a[0]>i or a[1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
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

## Poblaciones condicionadas

def population_conditional_array_alpha(array,aa1,aa2):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Obtiene de una array con los valores de ángulos diedros para un aminoácido una matriz de tamaño 360x360 con la población de cada uno de los puntos condicionada a que 
    el aminoácido vecino (aa2) tenga una conformación alfa ((-80,-20)+-50).
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de un aminóacidos de la cadena.
        Salida: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente, si el punto del aminoácido vecino presenta una estructura alfa.
    """
    array1=copy.deepcopy(array)
    array1[:,2*aa1]=array1[:,2*aa1]+180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los ángulos correspondientes
        k=0
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):           # Buscamos los índices que se identifican con los ángulos psi y phi en la fila a 
                if (a[2*aa2]>= -130 and a[2*aa2]<= -30) and (a[2*aa2+1]>= -70 and a[2*aa2+1]<=30):  # Incluimos la restricción de conformación alfa en el aminoácido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1

    return empty_array


def population_conditional_array_ppii(array,aa1,aa2):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Obtiene de una matriz con los valores de ángulos diedros para un aminoácido una matriz de tamaño 360x360 con la población de cada uno de los puntos condicionada a que 
    el aminoácido vecino (aa2) tenga una conformación ppII ((-60,145)+-30).
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de un aminóacidos de la cadena.
        Salida: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente, si el punto del aminoácido vecino presenta una estructura ppII.
    """
    array1=copy.deepcopy(array)  # Para no modificar la matriz original asignamos una copia de la matriz en otra variable
    array1[:,2*aa1]=array1[:,2*aa1]+180             
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los ángulos correspondientes
        k=0
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):           # Buscamos los índices que se identifican con los ángulos psi y phi en la fila a 
                if (a[2*aa2]>= -90 and a[2*aa2]<= -30) and (a[2*aa2+1]>= 115 and a[2*aa2+1]<=175):  # Incluimos la restricción de conformación ppii en el aminoácido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:                               # Si los índices todavía no se asocian a los valores de los ángulos diedros seguimos recorriendo i y 
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1
    return empty_array


def population_conditional_array_beta(array,aa1,aa2):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Obtiene de una matriz con los valores de ángulos diedros para un aminoácido una matriz de tamaño 360x360 con la población de cada uno de los puntos condicionada a que 
    el aminoácido vecino (aa2) tenga una conformación beta ((-130,140)+-40).
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de un aminóacidos de la cadena.
        Salida: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente, si el punto del aminoácido vecino presenta una estructura beta.
    """
    array1=copy.deepcopy(array)  # Para no modificar la matriz original asignamos una copia de la matriz en otra variable
    array1[:,2*aa1]=array1[:,2*aa1]+180            # La primera fila de la matriz representa el ángulo -180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180     # La primera columna de la matriz representa el ángulo 180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los ángulos correspondientes
        k=0
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):           # Buscamos los índices que se identifican con los ángulos psi y phi en la fila a 
                if (a[2*aa2]>= -170 and a[2*aa2]<= -90) and (a[2*aa2+1]>= 120 and a[2*aa2+1]<=180):  # Incluimos la restricción de conformación beta en el aminoácido posterior o anterior
                    empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:                 # Si los índices todavía no se asocian a los valores de los ángulos diedros seguimos recorriendo i y j
                if a[2*aa1]>i:
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1
    return empty_array

def population_conditional_array_remaining(array,aa1,aa2):
    """
    (numpy.ndarray) -> (numpy.ndarray)
    Obtiene de una matriz con los valores de ángulos diedros para un aminoácido una matriz de tamaño 360x360 con la población de cada uno de los puntos condicionada a que 
    el aminoácido vecino (aa2) no presente ninguna de las conformaciones anteriormente descritas.
        Entrada: Una array formada por los datos de ángulos diedros (phi y psi) de un aminóacidos de la cadena.
        Salida: Una array de tamaño 360 x 360 en el que cada valor de la array representa la población de cada punto, que contienen la pareja de valores para
        ángulo phi y psi correspondiente, si el punto del aminoácido vecino presenta una estructura beta.
    """
    array1=copy.deepcopy(array)    # Para no modificar la matriz original asignamos una copia de la matriz en otra variable
    array1[:,2*aa1]=array1[:,2*aa1]+180             # La primera fila de la matriz representa el ángulo -180
    array1[:,2*aa1+1]=-1*array1[:,2*aa1+1]+180      # La primera columna de la matriz representa el ángulo 180
    
    empty_array=np.zeros((360,360),dtype=int)  # array 360x360 vacia
    for a in array1:
        i=0                                     # Angulo phi
        j=0                                     # Angulo psi                                     # Indica que se han encontrado los ángulos correspondientes
        k=0
        while (a[2*aa1]>i or a[2*aa1+1]>j or k==0):       # Recorremos el array hasta obtener los valores de ángulos psi y phi de cada fila
            if (a[2*aa1]<=i and a[2*aa1+1]<=j):           # Buscamos los índices que se identifican con los ángulos psi y phi en la fila a 
                if not(a[2*aa2]>= -130 and a[2*aa2]<= -30) or not(a[2*aa2+1]>= -70 and a[2*aa2+1]<=30):       # Excluimos conformación alpha
                    if not(a[2*aa2]>= -90 and a[2*aa2]<= -30) or not(a[2*aa2+1]>= 115 and a[2*aa2+1]<=175):  # Excluimos conformación ppii
                        if not(a[2*aa2]>= -170 and a[2*aa2]<= -90) or not(a[2*aa2+1]>= 120 and a[2*aa2+1]<=180):  # Excluimos conformación beta 
                            empty_array[i-1,j-1]=empty_array[i-1,j-1]+1
                k=1
            else:                                          # Si los índices todavía no se asocian a los valores de los ángulos diedros seguimos recorriendo i y j
                if a[2*aa1]>i:                          
                    i=i+1
                if a[2*aa1+1]>j:
                    j=j+1
    return empty_array

def alpha_ind_prob(array):
    """
    Calcula la probabilidad de que el aminoácido adopte una conformación ppII.
        Entrada: Un array con la población en cada uno de los valores phi y psi del aminoácido
        Salida: Probabilidad de adoptar probabilidad alfa en ese aminoácido
    """
    phi=-80+180             # Valores de los ángulos diedros correspondientes con la conformación alfa
    psi=20+180
    n=0
    for i in range(phi-50,phi+50):
        for j in range(psi-50,psi+50):
            n=n+array[i-1][j-1]
    return(n)

def ppii_ind_prob(array):
    """
    Calcula la probabilidad de que el aminoácido adopte una conformación ppII.
        Entrada: Un array con la población en cada uno de los valores phi y psi del aminoácido
        Salida: Probabilidad de adoptar probabilidad ppII en ese aminoácido
    """
    phi=-60+180             # Valores de los ángulos diedros correspondientes con la conformación ppII
    psi=-145+180
    n=0
    for i in range(phi-30,phi+30):
        for j in range(psi-30,psi+30):
            n=n+array[i-1][j-1]
    return(n)

def beta_ind_prob(array):
    """
    Calcula la probabilidad de que el aminoácido adopte una conformación beta.
        Entrada: Un array con la población en cada uno de los valores phi y psi del aminoácido
        Salida: Probabilidad de adoptar probabilidad beta en ese aminoácido
    """
    phi=-130+180            # Valores de los ángulos diedros correspondientes con la conformación beta
    psi=-140+180    
    n=0
    for i in range(phi-40,phi+40):
        for j in range(psi-40,psi+40):
            n=n+array[i-1][j-1]
    return(n)

def remaining_ind_prob(array):
    """
    Calcula la probabilidad de que el aminoácido no adopte ninguna de las conformaciones anteriores.
        Entrada: Un array con la población en cada uno de los valores phi y psi del aminoácido
        Salida: Probabilidad de adoptar probabilidad "remaining" en ese aminoácido
    """
    n=1-(alpha_ind_prob(array)+ppii_ind_prob(array)+beta_ind_prob(array))       # Probabilidad de no encontrar ninguna de las conformaciones anteriores
    return(n)

def main():
    print("Módulo con las funciones utilizadas para el transcurso del proyecto de Trabajo de Fin de Máster")

if __name__ == "__main__":
    main()