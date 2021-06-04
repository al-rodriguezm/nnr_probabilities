#!/bin/bash

#################################
# $1 Directorio de los ficheros #
# $2 Nombre de los ficheros     #     
#################################


# Eliminamos los directorios de trabajos previos
rm -rf arrays
rm -rf population
rm -rf plots
rm -rf probabilities

# Creamos una carpeta auxiliar donde realizar los distintos procesos
mkdir aux

cp  $1/tray0*/pent* aux

gunzip aux/$2*
cat aux/* | sed '/^ *$/d' > rama.dat

rm -rf aux

# Creamos las carpetas necesarias para el análisis del péptido actual

mkdir arrays 
mkdir population
mkdir plots
mkdir probabilities



